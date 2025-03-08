#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>

#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

int vector_len;       // length of vector we will sum

int num_threads;            // number of threads we will use
int num_processors;         // number of processors we will use
pthread_t* tid;             // array of thread ids

typedef struct {
    long i;
    signal* sig;
    double bandwidth;
    int filter_order;
    double* filter_coeffs;
    double* band_power;
  } WorkerArgs;

void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands\n");
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}

void* worker(void* arg) {
    WorkerArgs* worker_args = (WorkerArgs*)arg;
    long myid     = worker_args->i;
    int blocksize = vector_len / num_threads; // note: floor
  
    // put ourselves on the desired processor
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(myid % num_processors, &set);
    if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
      perror("Can't setaffinity"); // hopefully doesn't fail
      exit(-1);
    }
  
    // This figures out the chunk of the vector I should
    // work on based on my id
    int mystart = myid * blocksize;
    int myend   = 0;
    if (myid == (num_threads - 1)) { // last thread
      // the last thread will take care of the leftover
      // elements of the vector, in case num_threads doesn't
      // divide vector_len
      // WARNING: this is a suboptimal solution. It means that the last thread
      // might do much more work than the other threads (up to almost double)
      // which will slow down the entire job. A better solution would split up
      // remainder work equally between threads...
      myend = vector_len;
    } else {
      myend = (myid + 1) * blocksize;
    }
    printf("Start and end of thread %ld: %d to %d\n", myid, mystart, myend);
    for (int j = mystart; j < myend; j++) {
        generate_band_pass(worker_args->sig->Fs,
            j * worker_args->bandwidth + 0.0001, // keep within limits
            (j + 1) * worker_args->bandwidth - 0.0001,
            worker_args->filter_order,
            worker_args->filter_coeffs);
        hamming_window(worker_args->filter_order,worker_args->filter_coeffs);

        // Convolve
        convolve_and_compute_power(worker_args->sig->num_samples,
            worker_args->sig->data,
            worker_args->filter_order,
            worker_args->filter_coeffs,
                    &(worker_args->band_power[j]));
    }
  
    // Done.  The master thread will sum up the partial sums
    pthread_exit(NULL);           // finish - no return value
  
}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub) {

  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  double filter_coeffs[filter_order + 1];
  double band_power[num_bands];

  tid = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);

  WorkerArgs all_worker_args[num_threads];

  for (int i = 0; i < num_threads; i++) {

    all_worker_args[i].i = i;
    all_worker_args[i].sig = sig;
    all_worker_args[i].bandwidth = bandwidth;
    all_worker_args[i].filter_order = filter_order;
    all_worker_args[i].filter_coeffs = malloc(sizeof(double) * (filter_order + 1));
    memcpy(all_worker_args[i].filter_coeffs, filter_coeffs, (filter_order + 1) * sizeof(double));
    all_worker_args[i].band_power = band_power;
    
    int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                    NULL, // use default attributes
                                    worker, // thread will begin in this function
                                    (void*)(&all_worker_args[i]) // we'll give it i as the argument
                                    );
    if (returncode != 0) {
      perror("Failed to start thread");
      exit(-1);
    }
  }

  for (int i = 0; i < num_threads; i++) {
    int returncode = pthread_join(tid[i], NULL);
    if (returncode != 0) {
      perror("join failed");
      exit(-1);
    }
    free(all_worker_args[i].filter_coeffs);
  }


  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
User time        %lf seconds\n\
System time      %lf seconds\n\
Page faults      %ld\n\
Page swaps       %ld\n\
Blocks of I/O    %ld\n\
Signals caught   %ld\n\
Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  num_threads = atoi(argv[6]); // number of threads
  num_processors = atoi(argv[7]); // numer of processors to use

  if (num_threads > num_bands) {
    num_threads = num_bands;
  }

  vector_len = num_bands;

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}

