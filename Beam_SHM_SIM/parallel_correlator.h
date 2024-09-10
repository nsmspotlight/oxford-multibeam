/* Program for offline processing -- Written by Jayanta Roy (Last modified on 1 June 2007)*/
#ifndef CORRELATOR_H
#define CORRELATOR_H
/* header file correlator.h */
#include "gmrt_newcorr.h"

#define DEBUG
#undef DEBUG

#define GATHER

#ifdef COMP_BEAM
#define BEAM_MODE1
#define BEAM_MODE2
#endif

// #define CHANNEL 2048
// #define FFTLEN (4*CHANNEL)
// #define NCHAN 4
#define CORRLEN (NCHAN * FFTLEN)
#define NUM_POL 32 // No of pols computation is handling
// #define NUM_ACQ 16     // No of acquation node
#define NUM_PROCESSES 48 // Total no of nodes take part
// #define NUM_ANT 8        //  No of antennas take part into correlation

#ifdef POLAR_MODE
#define NPOL 2            // Increase of polar term per baseline comapre to Intensity mode
#define NPOL_T (2 * NPOL) // Total no of polar term per beseline
#define NUM_CORR 32       // No of correlation node for a given pol
#else
#define NPOL 1
#define NPOL_T 1
#define NUM_CORR 16 // No of correlation node for a given pol
#endif

#define FRNG_STEP 0.25
#define NSTEP 2880
// #define BandWidth 33333333

#define NCORR (NPOL_T * NUM_ANT * NCHAN * (NUM_ANT * NCHAN + 1) / 2)
// #define M_PI 3.141592654
#define ACQ_LEN (32 * 1024 * 1024)
#define UNPACK_LEN (256 * 1024 * 1024)
#define MPI_BUF_CHUNK (ACQ_LEN / NUM_ACQ)
#define MPI_EXCESS_CHUNK (64 * 1024)
#define MPI_OVL_CHUNK (MPI_BUF_CHUNK + MPI_EXCESS_CHUNK)
#define ACQ_OVL (ACQ_LEN + MPI_EXCESS_CHUNK)
#define BEAM_SIZE (UNPACK_LEN / (NCHAN * NUM_ANT))
#define BEAM_SCALE 2 * (pow(10, 11))

#define SPLIT1 3
#define SPLIT2 4

#define SWEEPS_PER_DUMP 1
#define CORR_SIZE ((FFTLEN * NCORR) / 2) //%2% Same size as for one pol!

#define NTCHUNK 16
#define FFTBLOCK 64

//static float corrbuf[2 * CORR_SIZE]; //%2% CorrSize become double
//static float corr_sum[2 * CORR_SIZE];

//#define NSLW 1
//#define SPLITSEND NSLW

//#ifdef SPLITSEND
//static float corr_send[2 * CORR_SIZE];
//#endif

//static float rdelay_ti[NCHAN * NUM_ACQ];
//static float phase_ti[2][NUM_ACQ * NCHAN * FFTLEN / 2];
//static float dphase_t0[2][NUM_ACQ * NCHAN * FFTLEN / 2];
//static float rphase_ti[NUM_ACQ * NCHAN * FFTLEN / 2];
//static char pabuf[2][BEAM_SIZE];
//static char pa_voltage[2][BEAM_SIZE * NUM_ANT];

//static short int iabuf[2][(BEAM_SIZE)*NUM_ACQ]; // 16*beam_len*4
//static short int iabuf1[(BEAM_SIZE)*NUM_ACQ];   // 16*beam_len*4

//typedef float corrarray_t[NCORR][FFTLEN];
//typedef float datarray_t[NCHAN * NUM_ANT][FFTLEN];

//static corrarray_t corrampphas;
//static corrarray_t corrbufavg;
//static corrarray_t corrbuf_acc;

void corr_driver(signed char *buffer, int nint_corr, int nint_ia, int nint_pa, struct timeval timestamps, double *ch_freq, double BW, int w_buff, int fstop, int PAPols, int IAPols, int IASubPols, float step, short int *phase, int *iamask, short int *pamask, int ia_count, int pa_count, FILE *fp_data, FILE *fp_time, float *phastab, char *beam_mode1, char *beam_mode2, int mode);
//void correlate(short int obuf[FFTLEN / FFTBLOCK][UNPACK_LEN / (FFTLEN * NCHAN * NUM_ANT * NPOL)][NUM_ANT * NPOL][NCHAN][FFTBLOCK], float corrbuf[FFTLEN / 16][NCORR][2][4], short int outiabuf[BEAM_SIZE / 4], short int outpabuf[(BEAM_SIZE * NPOL) / 4], char pabuf[BEAM_SIZE], int iamask[NCHAN * NUM_ACQ], short int pamask[NCHAN * NUM_ACQ], int nint_ia, int nint_pa, int ia_count, int pa_count);
//void convcorr(float corrbuf[CORR_SIZE], float outcorrbuf[FFTLEN / 8][NCORR][2][4]);
void write_corr(float *corr_buf, struct tm *local_t, double time_ms, struct timeval timestamps, int iteration, FILE *fp_data, FILE *fp_time);
//void write_beam(short int *iabuf, int nint_ia);

/* corrdas.h methods oct4, 2010
float current_time(void);
int send_data(int integrated, double time_ms, float *corr_buf);
int corr_hdr(char *filename, CorrType *corr);
int corr_pre(char *filename, CorrType *corr);
void get_antenna(FILE *f, AntennaParType *antenna);
int get_sampler(FILE *f, CorrType *corr);
void get_freq_ch0_multi(char *gtac_code, SourceParType *source)
void get_freq_ch0(SourceParType *source);
double mjd_cal(struct tm *t);
double lmst(double mjd);
void calModel(double tm, CorrType *corr, ModelParType *par, int scans,ScanInfoType *scaninfo, AntennaFlagType *antflag);
void calModelPar(double tm, CorrType *corr, ModelParType *mpar, SourceParType *source, int *antmask);
void gsbe_ScanInfo(CorrType *corr, ModelInfoType *model, ScanInfoType *) ;  GSBC : Just getting info for each scan*/

// ================ Variables declration for acqpsr routinces.. ======================//
int BlocksRead;
// =================   function declration for acqpsr routinces.. ======================//
int initialize_vlt_shm();
// int initialize_IAPA(int pols, int curr_chans, int beam_integration, int output_size, int beamid);
// int initialize_IAPA(int pols, int curr_chans, int beam_integration, int output_size, int ddc, int decimation, int beamid, int nBMs);
int initialize_IAPA(int pols, int curr_chans, int beam_integration, int output_size, int beamid, int nBMs);
// int WriteBeam2SHM(char *beam_mode, int IAPols,int IASubPols,int PAPols, int ia_integ, int pa_integ, struct timeval timestamp, struct timeval timestamp_gps, double blk_nano, short int *iab, short int *col_pa, int iteration);
int WriteBeam2SHM(int nstokes, struct timeval timestamp, struct timeval timestamp_gps, double blk_nano, char *iab, int iteration, long beam_blocksize, int nBMs);
int WriteVltBeam2SHM(struct timeval *timestamp, char *vlt_beam_data, long beam_blocksize, int pol, double blk_nano, int iteration);
int process(int rec_num);
int set_seq_num(int rec_num, double pc_time);
void set_blk_times(double pc_time, int seq, int time_ok);
void init_blk_eqn(int rec_num);
int get_rec(int rec_num, int accept, struct timeval timestamp);
int set_rec_time(int rec_num, struct timeval timestamp);
void get_intr_times(void);
int set_blk_eqn(void);
void blk_log(void);
void seq_log(int rec_num);
int init_gps(void);
double test_gps(int ind2);
int set_gps(int i);
void shift(int ind, int step, int last, int max, double *val);
void ist_log(void);
// ====================== function declration for acqpsr routinces.. ==========================//

#endif
