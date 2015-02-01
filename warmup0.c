#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define SQRT_PI 1.77245385090552
#define PI 3.14159265358979323846264338327950


void err_sys(const char* x) { 
    perror(x); 
    exit(1); 
}


void shift_data(gsl_complex_packed_array* data_ptr,int N_samples){
  double tmp;
  int i;
  for(i=0;i<N_samples/2;i++){
    tmp = REAL(*data_ptr,i);
    REAL(*data_ptr,i) = REAL(*data_ptr,i+N_samples/2);
    REAL(*data_ptr,i+N_samples/2) = tmp;
    tmp = IMAG(*data_ptr,i);
    IMAG(*data_ptr,i) = IMAG(*data_ptr,i+N_samples/2);
    IMAG(*data_ptr,i+N_samples/2) = tmp;
  }
  printf("INFO: shifted data\n");
  fflush(stdout);
}

void scale_data(gsl_complex_packed_array* data_ptr,int N_samples, double scale){
  int i;
  for(i=0;i<N_samples;i++){
    REAL(*data_ptr,i) *= scale;
    IMAG(*data_ptr,i) *= scale;
  }
  printf("INFO: scaled data by %f\n", scale);
  fflush(stdout);
}

void fprint_data(FILE* stream,gsl_complex_packed_array data,int N_samples, int freq_domain_flag){
  int i;
  double tf,dtf;                    // tf carries x-axis time or frequency
  if(freq_domain_flag){            // frequency domain
    tf =  -(double)N_samples/2.0; // Nyquist critical freq
    dtf = 1.0;                      // Maximum time is frequency step
  }else{                           // time domain
    tf = 0.0;                       // assumes time starts at zero
    dtf = 1.0/(N_samples);
  }
  for(i=0;i<N_samples;i++,tf+=dtf){
    if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n", tf, REAL(data,i), IMAG(data,i)))<0)
      err_sys("fprintf error");
  }
  if(freq_domain_flag) // freq always periodic, more debuggable this way
    if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n", tf, REAL(data,0), IMAG(data,0)))<0)
      err_sys("fprintf error");
}

void plot(gsl_complex_packed_array data, int N_samples, int freq_domain_flag){
  FILE * gnuplot_pipe;
  int i, gnuplot_lines = 2;
  char * commandsForGnuplot[] = {"set notitle\n",
            "plot '-' u 1:2 with lines notitle, '-' u 1:3 with lines notitle\n"};
  printf("INFO: invoking gnuplot\n");
  fflush(stdout);
  if((gnuplot_pipe = popen("gnuplot -persistent", "w")) == NULL)
    err_sys("popen in plot error");
  for (i=0;i<gnuplot_lines; i++){ //Send commands to gnuplot one by one.
    if((fprintf(gnuplot_pipe, "%s \n", commandsForGnuplot[i]))<0)
      err_sys("fprintf in plot error");
  }
  // once for real values
  fprint_data(gnuplot_pipe, data, N_samples,freq_domain_flag);
  fprintf(gnuplot_pipe, "e");
  // And once for complex values
  fprint_data(gnuplot_pipe, data, N_samples,freq_domain_flag);
  fprintf(gnuplot_pipe, "e");
  fflush(gnuplot_pipe);
  fclose(gnuplot_pipe);
}

void plot_2_curves(gsl_complex_packed_array data0, gsl_complex_packed_array data1, int N_samples, int freq_domain_flag){
  FILE * gnuplot_pipe;
  int i, gnuplot_lines = 2;
  char * commandsForGnuplot[] = {"set notitle\n",
            "plot '-' u 1:2 with lines title 'Data0 Real', "
                 "'-' u 1:3 with lines title 'Data0 Complex', "
                 "'-' u 1:2 with lines title 'Data1 Real', "
                 "'-' u 1:3 with lines title 'Data1 Complex'\n"};
  printf("INFO: invoking gnuplot, 2 arrays\n");
  fflush(stdout);
  if((gnuplot_pipe = popen("gnuplot -persistent", "w")) == NULL)
    err_sys("popen in plot error");
  for (i=0;i<gnuplot_lines; i++){ //Send commands to gnuplot one by one.
    if((fprintf(gnuplot_pipe, "%s \n", commandsForGnuplot[i]))<0)
      err_sys("fprintf in plot error");
  }
  // once for real values data0
  fprint_data(gnuplot_pipe, data0, N_samples,freq_domain_flag);
  fprintf(gnuplot_pipe, "e");
  // And once for complex values data0
  fprint_data(gnuplot_pipe, data0, N_samples,freq_domain_flag);
  fprintf(gnuplot_pipe, "e");
  // once for real values data1
  fprint_data(gnuplot_pipe, data1, N_samples,freq_domain_flag);
  fprintf(gnuplot_pipe, "e");
  // And once for complex values data1
  fprint_data(gnuplot_pipe, data1, N_samples,freq_domain_flag);
  fprintf(gnuplot_pipe, "e");
  fflush(gnuplot_pipe);
  fclose(gnuplot_pipe);
}

void save_data(char* filename, gsl_complex_packed_array data, int N_samples, int freq_domain_flag){
  FILE* f;
  if((f = fopen(filename, "w")) == NULL) err_sys("fopen in save_data error");
  fprint_data(f,data,N_samples, freq_domain_flag);
  fclose(f);
  printf("INFO: saved data in file %s\n",filename);
  fflush(stdout);
}

void fft(double ** data_ptr, int N, int *freq_domain_flag_ptr){
  printf("INFO: Transforms data\n");
  fflush(stdout);
  if(gsl_fft_complex_radix2_forward(*data_ptr,1,N) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_forward error");
  *freq_domain_flag_ptr = !(*freq_domain_flag_ptr);
}

void ifft(double ** data_ptr, int N, int *freq_domain_flag_ptr){
  fflush(stdout);
  if(gsl_fft_complex_radix2_inverse(*data_ptr,1,N) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_inverse error");
  *freq_domain_flag_ptr = !(*freq_domain_flag_ptr);
  printf("INFO: Inverse transformed data\n");
}

void print_horror_diffs(gsl_complex_packed_array data0, gsl_complex_packed_array data1, int N, double eps){
  int i,flag=0;
  for(i=0;i<N;i++){
    if(fabs(REAL(data0,i) - REAL(data1,i)) > eps){
      printf("Horror diff: % 9.7f Real\n", REAL(data0,i)-REAL(data1,i));
      flag=1;
    }
    if(fabs(IMAG(data0,i) - IMAG(data1,i)) > eps){
      printf("Horror diff: % 9.7f Imag\n", IMAG(data0,i)-IMAG(data1,i));
      flag=1;
    }
  }
  if(!flag) printf("INFO: searched and found no horror diffs\n");
}

int main(int argc, char *argv[]){
  gsl_complex_packed_array data = calloc(2048,sizeof(double));
  double compar[2048];
  double dt = 1.0/1024.0, sigma = 1.0/64.0, sipisig,dtdt,sigmasigma,pisigmapisigma;
  int i,N = 1024,flag=0;
  sipisig = 1.0/(SQRT_PI*sigma);
  sigmasigma = sigma*sigma;
  pisigmapisigma=sigmasigma*PI*PI;
  dtdt = dt*dt;

  // GSL puts t = 0 at element 0.
  // For something to be periodic we need to set values from
  // right and left simultanously
  for(i=0;i<N/2;i++){
    REAL(data,i) = sipisig*exp(-i*i*dtdt/sigmasigma);
    REAL(data,N-1-i) = sipisig*exp(-(i+1)*(i+1)*dtdt/sigmasigma);
  }

  // Initializing compar
  for(i=0;i<N;i++){
    REAL(compar,i) = exp(-pisigmapisigma*(i-N/2)*(i-N/2));
  }

  save_data("before_fft.data",data,N,flag);
  fft(&data,N,&flag);
  save_data("after_fft.data",data,N,flag);

  //shift_data(&data,N);
  //scale_data(&data,N,1.0/(double)N);
  //plot_2_curves(data,compar,N,1);
  free(data);
  return 0;
}

