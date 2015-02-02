#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_fft_complex.h>
#include <pthread.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define MAX_CHARS_PER_LINE 100

struct fft_collection {
  gsl_complex_packed_array data; // Always allocate for log2_ceil(2*N) doubles
  unsigned char flag;   // Are we in freq domain and is data shifted?
  int N;                // actual found elements in data file. GSL want log2_ceil of this
  pthread_mutex_t lock; // Wrap for-loops over data elements in lock/unlock
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

void fft(struct fft_collection* coll){
  pthread_mutex_lock(&(coll->lock));
  if(gsl_fft_complex_radix2_forward(coll->data,1,log2_ceil(coll->N)) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_forward error");
  coll->flag = ~(coll->flag); // FFT will both change domain and shift data
  pthread_mutex_unlock(&(coll->lock));
}

void ifft(struct fft_collection* coll){
  pthread_mutex_lock(&(coll->lock));
  if(gsl_fft_complex_radix2_inverse(coll->data,1,log2_ceil(coll->N)) != GSL_SUCCESS)
    err_sys("gsl_fft_complex_radix2_inverse error");
  coll->flag = ~(coll->flag); // FFT will both change domain and shift data
  pthread_mutex_unlock(&(coll->lock));
}

/* Applies a gaussian filter to data in coll.
 * Assumes steps between frequencies is 1 (max time in sample is 1.0-dt). 
 */
void gaussian_filter(struct fft_collection* coll, double fc, double bandwidth){
  int i;
  double ib2 = 1.0/(bandwidth*bandwidth), scale;

  pthread_mutex_lock(&(coll->lock));
  scale = exp(-0.5*ib2*(0 - fc)*(0 - fc));
  REAL(coll->data,0) *= scale;
  IMAG(coll->data,0) *= scale;
  for(i=1;i<coll->N/2;i++){
    scale = exp(-0.5*ib2*(i - fc)*(i - fc));
    REAL(coll->data,i)         *= scale;
    REAL(coll->data,coll->N-i) *= scale;
    IMAG(coll->data,i)         *= scale;
    IMAG(coll->data,coll->N-i) *= scale;
  }
  scale = exp(-0.5*ib2*(i - fc)*(i - fc));
  REAL(coll->data,coll->N/2) *= scale;
  IMAG(coll->data,coll->N/2) *= scale;
  pthread_mutex_unlock(&(coll->lock));
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

/* All knowledge about the signal encapsulated here.
 * The found string is put into the_string.
 * Blindly trust the_string is long enough. */
void decode(struct fft_collection* coll, char *the_string){
  int bi, wave_packets = 128; // Magic number from manual inspection of signal
  int samp, samples_per_bit = coll->N/wave_packets;
  int by, bytes_in_data = wave_packets/8;
  int samples_per_byte = samples_per_bit*8;
  double max = 0.0, limit, sensor, val;
  char acc;

  // go through half the data to find the max
  for(bi=0;bi<wave_packets/2;bi++){
    for(samp=0;samp<samples_per_bit;samp++){
      if(fabs(REAL(coll->data, (samples_per_bit*bi)+samp)) > max)
        max = fabs(REAL(coll->data, (samples_per_bit*bi)+samp));
    }
  }
  limit = 0.8*max; // Magic number from... Just guessing

  for(by=0;by<bytes_in_data;by++){ // Go through each byte
    acc = 0x00;                    // Accumulate found 1-bits into this byte
    for(bi=0;bi<8;bi++){           // Go through each bit
      sensor = 0.0;
      for(samp=0;samp<samples_per_bit;samp++){ // Find max sampled value in this bit
        val = fabs(REAL(coll->data, (samples_per_byte*by)+(samples_per_bit*bi)+samp));
        if(val > sensor)
          sensor = val;
      }
      if(sensor > limit)           // Then we have found a 1-bit. Add it into acc
        acc |= 1 << (7-bi);
    }
    the_string[by] = acc;          // Put constructed ASCII-character into the_string
  }
}

int main(int argc, char *argv[]){
  FILE *f;
  struct fft_collection coll = {NULL, 0x00, 0, PTHREAD_MUTEX_INITIALIZER}; 
  int i;
  double fc = 1024.0, bandwidth = 256.0;
  char* the_string;

  if(argc != 2){
    printf("Usage: %s filename\n", argv[0]);
    return 1;
  }

  // See if file is any good...
  if((coll.N = count_samples(argv[1])) <= 1){
    printf("Found only %d data points in %s\n", coll.N, argv[1]);
    return 1;
  }

  // Ok, initialize rest of fft_collection struct
  if((coll.data = calloc(log2_ceil(coll.N<<1),sizeof(double))) == NULL)
    err_sys("calloc error");

  // Read data from file. Allcate memory for real and imaginary parts
  if ((f = fopen(argv[1], "r")) == NULL) err_sys("fopen error");
  for(i=0;i<coll.N;i++){ // IMAG(data,i) = 0.0; set by calloc
    fscanf(f,"%lf",&(REAL(coll.data,i)));
  }
  fclose(f);

  // Manipulate data
  fft(&coll);
  gaussian_filter(&coll, fc, bandwidth);
  ifft(&coll);

  // Find string in data
  if((the_string = calloc(log2_ceil(coll.N)/8+1,sizeof(char))) == NULL)
    err_sys("calloc error");
  decode(&coll, the_string);
  printf("%s\n", the_string);
  free(the_string);

  free(coll.data);
  return 0;
}
