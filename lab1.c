#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define MAX_CHARS_PER_LINE 100

void err_sys(const char* x) { 
    perror(x); 
    exit(1); 
}

/* Rounds a up to the nearest power of 2
 * Works for 1 < a < 1 073 741 825 only.
 * bsr instruction reason for making this */
int log2_ceil(int a){
  register int val __asm__("edi");
  __asm__("decl %edi;bsr %edi,%edi;incl %edi;");
  return 1 << val;
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

/* Encapsulates knowledge about discrete x-axis values corresponding to data
 * Assumes samples evenly distributed in time
 * makes time units so that t_k \in [ 0.0, 1.0-(1/N)] 
 *                          f_k \in [-N/2, N/2]
 * 7 decimal points is precision of input data
 * Create decimation-in-frequency wrapper of this that re-sorts data if
 * needed...
 * Scale like data/sqrt(N_samples) if in freq domain and data should fit
 * in same plot as time-domain data
 * */
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

void print_data(gsl_complex_packed_array data,int N_samples,int freq_domain_flag){
  fprint_data(stdout,data,N_samples,freq_domain_flag);
}


/* will create or supersede file. */
void save_data(char* filename, gsl_complex_packed_array data, int N_samples, int freq_domain_flag){
  FILE* f;
  if((f = fopen(filename, "w")) == NULL) err_sys("fopen in save_data error");
  fprint_data(f,data,N_samples, freq_domain_flag);
  fclose(f);
  printf("INFO: saved data in file %s\n",filename);
  fflush(stdout);
}

/* Won't store data in a file, streams directly to gnuplot */
void plot_data(gsl_complex_packed_array data, int N_samples, int freq_domain_flag){
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

/* Won't store data in a file, streams directly to gnuplot */
void plot_2_data(gsl_complex_packed_array data0, gsl_complex_packed_array data1, int N_samples, int freq_domain_flag){
  FILE * gnuplot_pipe;
  int i, gnuplot_lines = 2;
  char * commands_for_gnuplot[] = {"set notitle\n",
            "plot '-' u 1:2 with lines title 'Data0 Real', "
                 "'-' u 1:3 with lines title 'Data0 Complex', "
                 "'-' u 1:2 with lines title 'Data1 Real', "
                 "'-' u 1:3 with lines title 'Data1 Complex'\n"};
  printf("INFO: invoking gnuplot, 2 arrays\n");
  fflush(stdout);
  if((gnuplot_pipe = popen("gnuplot -persistent", "w")) == NULL)
    err_sys("popen in plot error");
  for (i=0;i<gnuplot_lines; i++){ //Send commands to gnuplot one by one.
    if((fprintf(gnuplot_pipe, "%s \n", commands_for_gnuplot[i]))<0)
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

int count_samples(char *filename){
  FILE *f;
  int N_samples;
  char buf[MAX_CHARS_PER_LINE];
  double dummy;
  if ((f = fopen(filename, "r")) == NULL)
    err_sys("fopen in count samples error");
  N_samples = 0;
  while(fgets(buf, MAX_CHARS_PER_LINE, f) != NULL &&
        sscanf(buf,"%lf",&dummy)>0 &&
        buf[0] != '\n')
    N_samples++;
  fclose(f);
  return N_samples;
}

void fft(double ** data_ptr, int N, int *freq_domain_flag_ptr){
  if(gsl_fft_complex_radix2_forward(*data_ptr,1,N) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_forward error");
  *freq_domain_flag_ptr = !(*freq_domain_flag_ptr);
  printf("INFO: Transformed data\n");
  fflush(stdout);
}

void ifft(double ** data_ptr, int N, int *freq_domain_flag_ptr){
  if(gsl_fft_complex_radix2_inverse(*data_ptr,1,N) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_inverse error");
  *freq_domain_flag_ptr = !(*freq_domain_flag_ptr);
  printf("INFO: Inverse transformed data\n");
  fflush(stdout);
}

int main(int argc, char *argv[]){
  FILE *f;
  int N_samples, N_memory, i;
  gsl_complex_packed_array data;
  int freq_domain_flag = 0; // input data most often in time domain

  if(argc != 3){
    printf("Usage: %s filename stream|save|streamsave\n", argv[0]);
    return 1;
  }

  // See if file is any good...
  if((N_samples = count_samples(argv[1])) <= 0){
    printf("Found no data points in %s\n", argv[1]);
    return 1;
  }
  printf("INFO: %s found %d data points in %s\n", argv[0], N_samples, argv[1]);
  fflush(stdout);
  N_memory = log2_ceil(N_samples<<1);

  // Read data from file. Allcate memory for real and imaginary parts
  if ((f = fopen(argv[1], "r")) == NULL) err_sys("fopen error");
  if((data = calloc(N_memory,sizeof(double))) == NULL)
    err_sys("calloc error");
  for(i=0;i<N_samples;i++){ // IMAG(data,i) = 0.0; set by calloc
    fscanf(f,"%lf",&(REAL(data,i)));
  }

  // Manipulate data
  fft(&data,&freq_domain_flag,N_memory>>1);
  //ifft(&data,&freq_domain_flag,N_memory>>1);

  // What to do with manipulated data
  if(argv[2][1] == 'a')       // option was save
    save_data("gnuplot_temp.data", data, N_samples,freq_domain_flag);
  else if(argv[2][6] == '\0') // option was stream
    plot(data, N_samples,freq_domain_flag);
  else{                       // option was streamsave
    save_data("gnuplot_temp.data", data, N_samples,freq_domain_flag);
    plot(data, N_samples,freq_domain_flag);
  }

  fclose(f);
  free(data);
  return 0;
}
