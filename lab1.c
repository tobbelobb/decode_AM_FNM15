#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fft_complex.h>
#include <pthread.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
// Rightmost (least significant) bit is 0 (1) if we are in time (freq) domain
// Leftmost (second) bit is 0 (1) if we have un-shifted (shifted) data
// default state: 0x00 // 0b00000000
#define FREQ_DOMAIN 0x01 // 0b00000001
#define SHIFTED     0x02 // 0b00000010

// Usage of macros:
// TEST_STATE(coll->flag, FREQ_DOMAIN)   // zero if we're in time domain
// TEST_STATE(coll->flag, SHIFTED)       // zero if not shifted
// TOGGLE_STATE(coll->flag, FREQ_DOMAIN) // we changed domain (did Fourier Tr)
// TOGGLE_STATE(coll->flag, SHIFTED)     // we did Fourier Tr or shift_coll
#define   TEST_STATE(bbm, x) ((bbm) &  (x))
#define    SET_STATE(bbm, x) ((bbm) |= (x))
#define  CLEAR_STATE(bbm, x) ((bbm) &= ~(x))
#define TOGGLE_STATE(bbm, x) ((bbm) ^= (x))

#define SQRT_PI 1.77245385090552
#define PI 3.14159265358979323846264338327950
#define MAX_CHARS_PER_LINE 100

struct fft_collection {
  gsl_complex_packed_array data; // Always allocate for log2_ceil(2*N) doubles
  unsigned char flag;   // Are we in freq domain and are data shifted?
  int N;                // actual found elements in data file. GSL want log2_ceil of this
  pthread_mutex_t lock; // Always locked when data is read/changed
};

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

void shift_coll(struct fft_collection* coll){
  double tmp;
  int i;
  pthread_mutex_lock(&(coll->lock));
  for(i=0;i<(coll->N)/2;i++){
    tmp = REAL(coll->data,i);
    REAL(coll->data,i) = REAL(coll->data,i+(coll->N)/2);
    REAL(coll->data,i+(coll->N)/2) = tmp;
    tmp = IMAG(coll->data,i);
    IMAG(coll->data,i) = IMAG(coll->data,i+(coll->N)/2);
    IMAG(coll->data,i+(coll->N)/2) = tmp;
  }
  TOGGLE_STATE(coll->flag, SHIFTED);
  pthread_mutex_unlock(&(coll->lock));
  printf("INFO: shifted data\n");
  fflush(stdout);
}

void scale_coll(struct fft_collection* coll, double scale){
  int i;
  pthread_mutex_lock(&(coll->lock));
  for(i=0;i<coll->N;i++){
    REAL(coll->data,i) *= scale;
    IMAG(coll->data,i) *= scale;
  }
  pthread_mutex_unlock(&(coll->lock));
  printf("INFO: scaled data by %f\n", scale);
  fflush(stdout);
}

void square_coll(struct fft_collection* coll){
  int i;
  pthread_mutex_lock(&(coll->lock));
  for(i=0;i<coll->N;i++){
    REAL(coll->data,i) *= REAL(coll->data,i);
    IMAG(coll->data,i) *= IMAG(coll->data,i);
  }
  pthread_mutex_unlock(&(coll->lock));
  printf("INFO: squared data\n");
  fflush(stdout);
}

/* Encapsulates knowledge about discrete x-axis values corresponding to data
 * Assumes samples evenly distributed in time
 * makes time units so that t_k \in [ 0.0, 1.0-(1/N)] 
 *                          f_k \in [-N/2, N/2]
 * 7 decimal points is precision of input data
 * Create decimation-in-frequency wrapper of this that re-sorts data if
 * needed...
 * Scale like data/sqrt(N) if in freq domain and data should fit
 * in same plot as time-domain data.
 * If coll->flag is set then data should always be shifted,
 * this could have been done in fprint_coll, but then printing
 * will start before shifting finished!!
 * Could be mitigated if each coll got its own mutex...
 * Lock should probably be called before calling this in most cases
 * */
void fprint_coll(FILE* stream, struct fft_collection* coll){
  int i;
  double tf,dtf; // carries x-axis value through loop
  // want plot-data to be unshifted
  if(TEST_STATE(coll->flag, SHIFTED)) shift_coll(coll);
  // Want differend x-axis for freq-domain data and time-domain data
  if(TEST_STATE(coll->flag, FREQ_DOMAIN)){ // frequency domain
    tf =  -(double)coll->N/2.0;            // Nyquist critical freq
    dtf = 1.0;                             // Maximum time is frequency step
  }else{                                   // time domain
    tf = 0.0;                              // assumes time starts at zero
    dtf = 1.0/(double)coll->N;
  }
  pthread_mutex_lock(&(coll->lock));
  for(i=0;i<coll->N;i++,tf+=dtf){
    if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n",
                        tf,    REAL(coll->data,i), IMAG(coll->data,i)))<0)
      err_sys("fprintf error");
  }
  if(TEST_STATE(coll->flag,FREQ_DOMAIN)) // freq always periodic, more debuggable this way
    if((fprintf(stream,"% 9.7e % 9.7e % 9.7e\n",
            tf,    REAL(coll->data,0), IMAG(coll->data,0)))<0)
      err_sys("fprintf error");
  pthread_mutex_unlock(&(coll->lock));
}

void print_coll(struct fft_collection* coll){
  fprint_coll(stdout,coll);
}


void copy_coll_data(struct fft_collection* source, struct fft_collection* dest){
  int i;
  if(source->N == dest->N){
    pthread_mutex_lock(&(dest->lock));
    pthread_mutex_lock(&(source->lock));
    for(i=0;i<source->N;i++){
      REAL(dest->data, i) = REAL(source->data, i);
      IMAG(dest->data, i) = IMAG(source->data, i);
    }
    pthread_mutex_unlock(&(dest->lock));
    pthread_mutex_unlock(&(source->lock));
  }else{
    printf("INFO: Didn't copy data. Different coll.N\n");
  }
  printf("INFO: Copied data\n");
}

/* will create or supersede file. */
void save_coll(struct fft_collection* coll, const char* filename){
  FILE* f;
  if((f = fopen(filename, "w")) == NULL) err_sys("fopen in save_coll error");
  fprint_coll(f,coll);
  fclose(f);
  printf("INFO: saved data in file %s\n",filename);
  fflush(stdout);
}

/* Won't store data in a file unless non NULL filename is given
 * */
void plot_coll(struct fft_collection* coll, const char* filename){
  FILE * gnuplot_pipe;
  int i, gnuplot_lines = 2;
  char * commandsForGnuplot[] = {"set notitle\n",
            "plot '-' using 1:2 with lines title 'Real',"
                 "'-' using 1:3 with lines title 'Imag'\n"};
  // copy coll                       make space for own data           copy flag   copy N   create own mutex.
  struct fft_collection coll_copy = {calloc(2*coll->N,sizeof(double)), coll->flag, coll->N, PTHREAD_MUTEX_INITIALIZER};
  // make a copy of the data, so coll can go on with its life in another function
  copy_coll_data(coll,&coll_copy); 

  printf("INFO: invoking gnuplot\n");
  fflush(stdout);
  if((gnuplot_pipe = popen("gnuplot -persistent", "w")) == NULL)
    err_sys("popen in plot error");
  for (i=0;i<gnuplot_lines; i++){ //Send commands to gnuplot one by one.
    if((fprintf(gnuplot_pipe, "%s \n", commandsForGnuplot[i]))<0)
      err_sys("fprintf in plot error");
  }
  // once for real values (gnuplot won't store streamed data)
  fprint_coll(gnuplot_pipe, &coll_copy);
  fprintf(gnuplot_pipe, "e\n");
  fflush(gnuplot_pipe);
  // And once for complex values (since gnuplot won't store streamed data)
  fprint_coll(gnuplot_pipe, &coll_copy);
  fprintf(gnuplot_pipe, "e\n");
  fflush(gnuplot_pipe);
  fclose(gnuplot_pipe);
  // lastly, if we wanted to save, save the copy
  if(filename != NULL) save_coll(&coll_copy, filename);
  free(coll_copy.data);
}

void fft(struct fft_collection* coll){
  pthread_mutex_lock(&(coll->lock));
  if(gsl_fft_complex_radix2_forward(coll->data,1,log2_ceil(coll->N)) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_forward error");
  coll->flag = ~(coll->flag); // FFT will both change domain and shift data
  pthread_mutex_unlock(&(coll->lock));
  printf("INFO: Transformed data\n");
  fflush(stdout);
}

void ifft(struct fft_collection* coll){
  pthread_mutex_lock(&(coll->lock));
  if(gsl_fft_complex_radix2_inverse(coll->data,1,log2_ceil(coll->N)) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_inverse error");
  coll->flag = ~(coll->flag); // FFT will both change domain and shift data
  pthread_mutex_unlock(&(coll->lock));
  printf("INFO: Inverse transformed data\n");
  fflush(stdout);
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

int main(int argc, char *argv[]){
  FILE *f;
  struct fft_collection coll = {NULL, 0x00, 0, PTHREAD_MUTEX_INITIALIZER}; 
  int i;

  if(argc != 3){
    printf("Usage: %s filename [p][s]|letter\n", argv[0]);
    return 1;
  }

  // See if file is any good...
  if((coll.N = count_samples(argv[1])) <= 1){
    printf("Found only %d data points in %s\n", coll.N, argv[1]);
    return 1;
  }
  printf("INFO: %s found %d data points in %s\n", argv[0], coll.N, argv[1]);
  fflush(stdout);

  // Ok, initialize rest of fft_collection struct
  if((coll.data = calloc(log2_ceil(coll.N<<1),sizeof(double))) == NULL)
    err_sys("calloc error");

  // Read data from file. Allcate memory for real and imaginary parts
  if ((f = fopen(argv[1], "r")) == NULL) err_sys("fopen error");
  for(i=0;i<coll.N;i++){ // IMAG(data,i) = 0.0; set by calloc
    fscanf(f,"%lf",&(REAL(coll.data,i)));
  }

  // Manipulate data
  fft(&coll);
  //ifft(&coll);

  // What to do with manipulated data
  if((argv[2][0] == 's' && argv[2][1] == 'p') ||
     (argv[2][0] == 'p' && argv[2][1] == 's')) // plot and save 
    plot_coll(&coll,"plotted_and_saved.data");
  else if(argv[2][0] == 's') // save
    save_coll(&coll, "saved.data");
  else // plot
    plot_coll(&coll, NULL);

  fclose(f);
  free(coll.data);
  return 0;
}
