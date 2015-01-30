#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define SQRT_PI 1.77245385090552

void err_sys(const char* x) { 
    perror(x); 
    exit(1); 
}

void fprint_data(FILE* stream,gsl_complex_packed_array data,int N_samples, int freq_domain_flag){
  int i;
  double tf,dt,dtf,scale;
  gsl_complex_packed_array ordered_data;
  if(freq_domain_flag){
    ordered_data = calloc((N_samples<<1),sizeof(double));
    scale = 1.0/sqrt(N_samples);
    for(i=0;i<N_samples/2;i++){
      REAL(ordered_data,i)             = scale*fabs(REAL(data,i+N_samples/2));
      IMAG(ordered_data,i)             = scale*fabs(IMAG(data,i+N_samples/2));
      REAL(ordered_data,i+N_samples/2) = scale*fabs(REAL(data,i));
      IMAG(ordered_data,i+N_samples/2) = scale*fabs(IMAG(data,i));
    }
    data = ordered_data;
  }

  dt = 1.0/((double)N_samples);
  if(freq_domain_flag){  // we are in frequency domain
    //tf =  -1.0/(2.0*dt); // Nyquist critical freq
    tf =  -(double)N_samples/(2.0); // Nyquist critical freq
    dtf = 1.0;
  }else{                 // we are in time domain
    tf = -0.5;
    dtf = dt;
  }
  for(i=0;i<N_samples;i++,tf+=dtf){
    if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n", tf, REAL(data,i), IMAG(data,i)))<0)
      err_sys("fprintf error");
  }
  if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n", tf, REAL(data,0), IMAG(data,0)))<0)
    err_sys("fprintf error");
  if(freq_domain_flag){
    free(ordered_data);
  }
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

void save_data(char* filename, gsl_complex_packed_array data, int N_samples, int freq_domain_flag){
  FILE* f;
  if((f = fopen(filename, "w")) == NULL) err_sys("fopen in save_data error");
  fprint_data(f,data,N_samples, freq_domain_flag);
  fclose(f);
  printf("INFO: saved data in file %s\n",filename);
  fflush(stdout);
}

void fft(double ** data_ptr, int *freq_domain_flag_ptr, int N){
  printf("INFO: Transforms data\n");
  fflush(stdout);
  if(gsl_fft_complex_radix2_forward(*data_ptr,1,N) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_forward error");
  *freq_domain_flag_ptr = !(*freq_domain_flag_ptr);
}

void ifft(double ** data_ptr, int *freq_domain_flag_ptr, int N){
  printf("INFO: Inverse transforms data\n");
  fflush(stdout);
  if(gsl_fft_complex_radix2_inverse(*data_ptr,1,N) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_inverse error");
  *freq_domain_flag_ptr = !(*freq_domain_flag_ptr);
}

int main(int argc, char *argv[]){
  gsl_complex_packed_array data = calloc(2048,sizeof(double));
  double dt = 1.0/1024.0, sigma = 1.0/64.0, sipisig,dtdt,sigmasigma;
  int i,N = 1024,flag=0;
  sipisig = 1.0/(SQRT_PI*sigma);
  sigmasigma = sigma*sigma;
  dtdt = dt*dt;

  for(i=0;i<N;i++){
    REAL(data,i) = sipisig*exp(-(-((double)N)/2.0+(double)i)*(-((double)N)/2.0+(double)i)*dtdt/sigmasigma);
  }

  save_data("warmup0_orig.data",data,N,flag);
  fft(&data,&flag,N);
  save_data("warmup0_transf.data",data,N,flag);
  free(data);
  return 0;
}
