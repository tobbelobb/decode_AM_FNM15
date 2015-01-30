#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define MAX_CHARS_PER_LINE 100

/* Rounds a up to the nearest power of 2
 * Works for 1 < a < 1 073 741 825 only.
 * bsr instruction reason for making this */
int log2_ceil(int a){
  register int val __asm__("edi");
  __asm__("decl %edi;bsr %edi,%edi;incl %edi;");
  return 0x0001 << val;
}

void err_sys(const char* x) { 
    perror(x); 
    exit(1); 
}

/* Encapsulates knowledge about discrete x-axis values corresponding to data
 * Assumes samples evenly distributed in time
 * makes time units so that t_k \in [ 0.0, 1.0] 
 *                          f_k \in [-1/(2\Delta), 1/(2\Delta)]
 * 7 decimal points is precision of input data
 * Create decimation-in-frequency wrapper of this that re-sorts data if
 * needed...
 * Scale like data/sqrt(N_samples) if in freq domain and data should fit
 * in same plot as time-domain data
 * */
void fprint_data(FILE* stream,gsl_complex_packed_array data,int N_samples, int freq_domain_flag){
  int i;
  double tf,dt,dtf,scale;
  gsl_complex_packed_array ordered_data;
  if(freq_domain_flag){ // Then we have unordered data...
    ordered_data = calloc(N_samples<<1,sizeof(double));
    scale = 1.0/sqrt(N_samples);
    for(i=0;i<N_samples/2;i++){
      REAL(ordered_data,i)             = scale*fabs(REAL(data,i+N_samples/2));
      IMAG(ordered_data,i)             = scale*fabs(IMAG(data,i+N_samples/2));
      REAL(ordered_data,i+N_samples/2) = scale*fabs(REAL(data,i));
      IMAG(ordered_data,i+N_samples/2) = scale*fabs(IMAG(data,i));
    }
    data = ordered_data;
  }

  dt = 1.0/((double)N_samples-1.0);
  if(freq_domain_flag){  // we are in frequency domain
    tf =  -1.0/(2.0*dt); // Nyquist critical freq
    dtf = 1.0;
  }else{                 // we are in time domain
    tf = -0.5;
    dtf = dt;
  }
  for(i=0;i<N_samples;i++,tf+=dtf){
    if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n", tf, REAL(data,i), IMAG(data,i)))<0)
      err_sys("fprintf error");
  }
  if(freq_domain_flag) free(ordered_data);
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
void plot(gsl_complex_packed_array data, int N_samples, int freq_domain_flag){
  FILE * gnuplot_pipe;
  int i, gnuplot_lines = 2;
  char * commandsForGnuplot[] = {"set notitle\n",
            "plot '-' u 1:2 with linespoints pointtype 7 notitle\n"};
  printf("INFO: invoking gnuplot\n");
  fflush(stdout);
  if((gnuplot_pipe = popen("gnuplot -persistent", "w")) == NULL)
    err_sys("popen in plot error");
  for (i=0;i<gnuplot_lines; i++){ //Send commands to gnuplot one by one.
    if((fprintf(gnuplot_pipe, "%s \n", commandsForGnuplot[i]))<0)
      err_sys("fprintf in plot error");
  }
  fprint_data(gnuplot_pipe, data, N_samples,freq_domain_flag);
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
  FILE *f;
  int N_samples,N_memory,i,freq_domain_flag=0; // input data most often in time domain
  gsl_complex_packed_array data;

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
