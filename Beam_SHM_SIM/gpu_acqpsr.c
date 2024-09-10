#include <stdio.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "gmrt_newcorr.h"
#include "acqpsr.h"
#include "pcidev2.h"
#include "parallel_correlator.h"
#include "gpu_corr_beam.h"

#define MAX_BLOCKS 8
#define BeamHdrKey 1050

//=================== ACQPSR RELARTED ============================== //
const double GPSmin = 60;
//static int ShmId = 0, DefBufs = 5;
static int ShmId = 0;
// static int       ia_ShmId=0, pa_ShmId=0;
static char *Shmp = NULL;
static GlobalInfoType *Info = NULL;
static RackInfoType RackInfo;
static RecType *Rec = NULL;
static long *BlockTime = NULL;
// static int       Integ=1, PAPols=2, IASubPols=2, IAPols=2;
static BlockType BlockInfo;
long RefDayTime = 0, TZ[2], TZoff = 0;
double BLKtime;
static double GPStime, GetIntr = 0;
static double BLKtimeEps = 2e-10, GPStimeEps = 6e-3;
// static double    BLKtimeEps=2e-3, GPStimeEps=6e-3;
static double BlktTimeEps = 1e-4, TickTime = 0.01;
double *ISTeqn, *BLKeqn, *GPSval, *BLKval;
static double GPSival[MaxGPS];
static TimeType GPSint[MaxGPS], BLKint[MaxBLK];
static int Gind = 0; // Indix to GPSint
static int GPS = 0, GPSind = 0, ISTeqnSet = 0;
int BLK = 1, BLKind = 0; // BLKeqnSet=0;
unsigned *GPSseq, *BLKseq;
static int *RecNum = 0, DataSeq = 0;
static int BeamFirstWrite = 1;
//static int PortInUse = 0, DataBlockSize = 0, BlocksPerRec = 0;
static int DataBlockSize = 0, BlocksPerRec = 0;
static int CurrRecs = MaxRecs;
static int MaxTimeErr = MaxRecs / 2, DataRecs = 0, NotStop = 1;
static long int CurrSamples = MaxSamples;
static long int CurrRecSize = RecSize;
static int count = 1;
//static int FD = 0, GPSfd;
static int GPSfd;
static char *GPSdev = "/dev/gps";
static int read_time = 1;
static struct timeval timestamp_gps_beg;
static struct timeval timestamp_pc_beg;
static unsigned int curBlock1 = 0, curRec1 = 0;
static double blk_nano_beg;
static int iteration_beg;
static double blk_nano_beg;

int file_mode_set = 0;

enum
{
  SeqTag,
  IntTag,
  BLKTag,
  ISTTag,
  FlagMask = 3
};
typedef struct
{
  unsigned short flag, rec_time[3]; // rec_time codes 48 bits of a double
  unsigned short dummy, rec_num;
  unsigned int seq_num;
} SeqLogType;
/*
typedef struct
{
  unsigned short flag, dummy;
  unsigned int int_num;
  double int_time;
} IntLogType;
*/
typedef struct
{
  unsigned short flag, pc_time[3]; // pc_time  codes 48 bits of a double
  unsigned blk_ref;
  float pc2blk;
} BLKlogType;
typedef struct
{
  unsigned short flag, pc_time[3]; // pc_time  codes 48 bits of a double
  unsigned ist_ref;
  float pc2ist;
} ISTlogType;
static unsigned short LogBuf[sizeof(SeqLogType) / sizeof(short)];
static SeqLogType *SeqLog = (SeqLogType *)LogBuf;
//static IntLogType *IntLog = (IntLogType *)LogBuf;
static BLKlogType *BLKlog = (BLKlogType *)LogBuf;
static ISTlogType *ISTlog = (ISTlogType *)LogBuf;
//static char *LogFile = "/tmp/acqbin.log";
//static FILE *Logfd = NULL;
static int LogEvents = 0; // Default was 1, set 0 to avoid kepping log, log file is not opened. To be opened to write before making it 1.
char errlog[128] = {0};

// static int IAstoke_integ=0, PAstoke_integ=0, stoke_integ=0;
static int stoke_integ = 0;
static PortInfoType port_info;

static DataHeader1 *dataHdr1;
static DataBuffer1 *dataBuf1;
// static double blk_nano=0.0;

//===================== ACQPSR RELATED =========================================//
int WriteVltBeam2SHM(struct timeval *timestamp, char *vlt_beam_data, long beam_blocksize, int pol, double blk_nano, int iteration)
{
  static struct timeval local_time;
  static double local_blk_nano, frac_sec;
  static struct tm *local_t;

	/*
  static FILE *acq_tm_itr = NULL;
  if (acq_tm_itr == NULL)
  {
    if ((acq_tm_itr = fopen("/tmp/tm_itr_acq_vlt.dat", "w")) == NULL)
    {
      perror("/tmp/tm_itr_acq_vlt.dat");
      return -1;
    }
  }
	*/

  curBlock1 = curRec1 % MaxDataBlocks;
  dataBuf1->blockSize = VltBeamBufferSize;
  // CAUTION See all flags are updated in (pol == 2) condition, so calling sequence
  // to this routine has to be pol1 first, and pol2 later. (from gmrt_corr_cuda.cu)
  if (pol == 1)
  {
    // Copying data for pol1 and common header.
    // fprintf(stderr, "#####POL1 CHAN 1000 = %d\n", *(vlt_beam_data+1000));
    // memcpy(dataBuf1->data_p1+curBlock1*(VltBeamBufferSize), vlt_beam_data, (VltBeamBufferSize)); // first pol data?
    memcpy(dataBuf1->data_p1 + curBlock1 * (VltBeamBufferSize), vlt_beam_data, (beam_blocksize)); // first pol data?
  }
  else if (pol == 2)
  {
    // Copying only data, header is updated when pol == 1 above. Need not do again the same.
    // memcpy(dataBuf1->data_p2+curBlock1*(VltBeamBufferSize), vlt_beam_data, (VltBeamBufferSize)); // Second pol data?
    memcpy(dataBuf1->data_p2 + curBlock1 * (VltBeamBufferSize), vlt_beam_data, (beam_blocksize)); // Second pol data?

    local_t = localtime((const time_t *)timestamp);
    frac_sec = timestamp->tv_usec;
    // frac_sec += blk_nano/1e6;
    fprintf(stderr, "VLTG TIME %02d:%02d:%02d.%06d%03d CUR_BLOCK = %d ITR = %d\n", local_t->tm_hour, local_t->tm_min, local_t->tm_sec, (int)frac_sec, (int)blk_nano, curBlock1, iteration);
    //fprintf(acq_tm_itr, "VLTG TIME %02d:%02d:%02d.%06d%03d CUR_BLOCK = %ld ITR = %d\n", local_t->tm_hour, local_t->tm_min, local_t->tm_sec, (int)frac_sec, (int)blk_nano, curBlock1, iteration);

    // fprintf(stderr, "VOLTAGE BEAM TIME (WT): Sec = %ld uSec = %06ld %03d\n", timestamp->tv_sec, timestamp->tv_usec,(int) blk_nano);
    // memcpy(&dataHdr1->timestamp+curBlock1*sizeof(timestamp), timestamp, sizeof(timestamp));
    // memcpy(&dataHdr1->timestamp_gps+curBlock1*sizeof(timestamp), timestamp, sizeof(timestamp)); // check if timestamp is right?
    memcpy(&(dataHdr1->timestamp[curBlock1]), timestamp, sizeof(struct timeval));
    memcpy(&(dataHdr1->blk_nano[curBlock1]), &blk_nano, sizeof(blk_nano));
    memcpy(&(dataHdr1->iteration[curBlock1]), &iteration, sizeof(iteration));

    // memcpy(&local_time, &dataHdr1->timestamp+curBlock1*sizeof(timestamp), sizeof(timestamp));
    memcpy(&local_time, &(dataHdr1->timestamp[curBlock1]), sizeof(struct timeval));
    memcpy(&local_blk_nano, &(dataHdr1->blk_nano[curBlock1]), sizeof(blk_nano));
    // fprintf(stderr, "VOLTAGE BEAM TIME (RD): Sec = %ld uSec = %06ld %03d\n", local_time.tv_sec, local_time.tv_usec,(int) local_blk_nano);

    dataBuf1->flag = BufReady;
    dataBuf1->curBlock = curBlock1;
    dataHdr1->refTime = timestamp->tv_sec + timestamp->tv_usec / 1e6;
    dataHdr1->dataTime = dataHdr1->refTime;
    dataBuf1->curRecord = curRec1++;
  }
  return 0;
}

int WriteBeam2SHM(int nstokes, struct timeval timestamp_pc, struct timeval timestamp_gps, double blk_nano, char *iab, int iteration, long beam_blocksize, int nBMs)
{
  //static double local_blk_nano, frac_sec;
  static double frac_sec;
  static struct tm *local_t;

/*
  static FILE *acq_tm_itr = NULL;
  if (acq_tm_itr == NULL)
  {
    if ((acq_tm_itr = fopen("/tmp/tm_itr_acq_BEAM.dat", "w")) == NULL)
    {
      perror("/tmp/tm_itr_acq_BEAM.dat");
      return -1;
    }
  }
*/

  //  struct tm *local_t;
  //  local_t = localtime(&timestamp.tv_sec);
  //  gps_ts_sec = local_t->tm_sec + local_t->tm_min*60 + local_t->tm_hour*3600;

  // fprintf(stdout,"BEAM SIZE >>>> = %d \n",beam_blocksize);
  // char *beam_mode = "IA";

  timestamp_gps_beg.tv_sec = timestamp_gps.tv_sec;
  timestamp_gps_beg.tv_usec = timestamp_gps.tv_usec;
  iteration_beg = iteration;
  blk_nano_beg = blk_nano;
  timestamp_pc_beg.tv_sec = timestamp_pc.tv_sec;
  timestamp_pc_beg.tv_usec = timestamp_pc.tv_usec;

  local_t = localtime((const time_t *)&timestamp_pc);
  frac_sec = timestamp_pc.tv_usec;
  // fprintf(stderr, "# %02d:%02d:%02d.%06d%03d ACQ ITERATION = %d RECNUM = %d\n", local_t->tm_hour, local_t->tm_min, local_t->tm_sec, (int)frac_sec, (int)blk_nano, iteration, *RecNum);
  fprintf(stderr, "# %02d:%02d:%02d.%06d%03d\n", local_t->tm_hour, local_t->tm_min, local_t->tm_sec, (int)frac_sec, (int)blk_nano);

  //    frac_sec += blk_nano/1e6;
  // SSK fprintf(stderr, "IAPA TIME %02d:%02d:%02d.%06d%03d ITR = %d\n",local_t->tm_hour,local_t->tm_min,local_t->tm_sec,(int) frac_sec, (int)blk_nano, iteration);
  // fprintf(stderr, "IAPA TIME %02d:%02d:%02d.%f\n",local_t->tm_hour,local_t->tm_min,local_t->tm_sec,frac_sec);
  // fprintf(acq_tm_itr, "IAPA TIME %02d:%02d:%02d.%06d%03d ITR = %d\n",local_t->tm_hour,local_t->tm_min,local_t->tm_sec,(int) frac_sec, (int)blk_nano, iteration);
  // fprintf(stderr, "IAPA TIME %02d:%02d:%02d.%06d%03d ITR = %d\n",local_t->tm_hour,local_t->tm_min,local_t->tm_sec,(int) frac_sec, (int)blk_nano, iteration);

  //    fprintf(stderr, "PA BEAM TIME (WT): Sec = %ld uSec = %06ld %03d\n", timestamp_pc.tv_sec, timestamp_pc.tv_usec,(int) blk_nano);

  /*
if(read_time == 1 ) {
  timestamp_gps_beg.tv_sec = timestamp_gps.tv_sec;
  timestamp_gps_beg.tv_usec = timestamp_gps.tv_usec;
  iteration_beg = iteration;
  blk_nano_beg=blk_nano;
  timestamp_pc_beg.tv_sec = timestamp_pc.tv_sec;
  timestamp_pc_beg.tv_usec = timestamp_pc.tv_usec;
  read_time = 0;
}
  */

  //     if(strcmp(beam_mode, "IA") == 0) stoke_integ = IAstoke_integ;
  //     else if(strcmp(beam_mode, "PA") == 0) stoke_integ = PAstoke_integ;

  stoke_integ += 1;
  if (stoke_integ == 1)
  {
    if (BeamFirstWrite)
    {
      //    set_extra_bufs(DefBufs); // Why do we need this?
      set_io(1);
      //    get_port_info(&port_info);
      port_info.block_size = IABeamBufferSize;
      // DataBlockSize = IABeamBufferSize/nint_ia;
      DataBlockSize = beam_blocksize;
      /*
          fprintf(stdout,"BEAM >>>> Data Block Size = %d expecting from GPU OUT\n",DataBlockSize);
          if (CurrRecSize % DataBlockSize != 0)
          {
      //      if(errlog['p']==1) fprintf(AcqErrLog,"CurrRecSize=%d not a multiple of block_size=%d\n", CurrRecSize, DataBlockSize);
              fprintf(stdout,"BEAM >>>> CurrRecSize=%d not a multiple of block_size=%d\n", CurrRecSize, DataBlockSize);
          } else {
           //if(strcmp(beam_mode, "IA") == 0) { BlocksPerRec = CurrRecSize/(IAPols*DataBlockSize);}
           //else if(strcmp(beam_mode, "PA") == 0) { BlocksPerRec = CurrRecSize/(PAPols*DataBlockSize);}
           BlocksPerRec = 1;
           //if(strcmp(beam_mode, "IA") == 0) { BlocksPerRec = 1;}
           //else if(strcmp(beam_mode, "PA") == 0) { BlocksPerRec = 1;}
      //    if(errlog['p']==1) fprintf(AcqErrLog,"Blocks Per Rec  = %d\n\n",BlocksPerRec);
          fprintf(stdout,"BEAM >>>> CURRRECSIZE = %d DATABLOCKSIZE = %d Blocks Per Rec  = %d\n\n",CurrRecSize, DataBlockSize, BlocksPerRec);
        }
      */
      BlocksPerRec = 1;
      //  fprintf (stderr, "## Beam_blocksize = %d    NBMs = %d RECNUM = %d ################## \n", beam_blocksize, nBMs, *RecNum);
      // memcpy(Rec[*RecNum].begp, iab, beam_blocksize*nBMs);
      memcpy(Rec[*RecNum].begp, iab, beam_blocksize);
      // fprintf(stdout,"crossed memcpy in iteration %d\n",*RecNum); // SSK COMMENTED
      //    perror("###########################   memcpy :");
      BlocksRead += 1;

      memcpy(&Rec[*RecNum].timestamp_gps, &timestamp_gps_beg, sizeof(struct timeval));
      memcpy(&Rec[*RecNum].timestamp_pc, &timestamp_pc_beg, sizeof(struct timeval));
      Rec[*RecNum].AcqSeqNo = iteration_beg;
      Rec[*RecNum].blk_nano = blk_nano_beg;
      Rec[*RecNum].pc_time = timestamp_pc.tv_sec + timestamp_pc.tv_usec / 1e6; // added for timestamp to come from GSB SSK/26nov2010
      get_rec(*RecNum, 1, timestamp_pc);
      init_blk_eqn(*RecNum);

      Rec->rec_flag = RackInfo.data_flag | GoodData;
      Rec->rec_time = ToIST(Rec->pc_time);
      Rec->rec_seq = DataSeq;

      RackInfo.blocks = 0;
      RackInfo.time_err = 0;

      DataSeq++;
      (*RecNum)++;
      Info->acq_flag &= ~Marked;

      BeamFirstWrite = 0;
      stoke_integ = 0;
      // if(strcmp(beam_mode, "IA") == 0) IAstoke_integ=0;
      // else if(strcmp(beam_mode, "PA") == 0) PAstoke_integ=0;
    }
    else
    {

      // BEAM     fprintf (stderr, "##########     RecNum = %d    ##################### \n", *RecNum);
      // fprintf (stderr, "############################### ****** second iteration \n");

      memcpy(Rec[*RecNum].begp, iab, beam_blocksize);

      // fprintf(stdout,"crossed memcpy in iteration %d\n",*RecNum); // SSK COMMENTED
      BlocksRead += 1;
      //  fprintf(stdout, "==================================================================================>>>>>>>>>>>>>>>> BEAM WRITTEN\n");

      //  check maximum time difference between two consecutive data block reads
      //    get_intr_times();

      //  Read only 'DataRecs' number of shared memory data blocks
      //  DataRecs=0 -> ifinite reading //next condition changed <= 24/11/2003 r.m.dabade
      if (NotStop && DataRecs != 0)
        NotStop = ++count != DataRecs;

      //  mark the current data block number for writing
      Info->acq_flag |= Marked;
      Info->mark_num++;

      //    Does much to set record time, insted a new RecStartTime is used to fill the time from
      //    GSB cluster (master time).

      memcpy(&Rec[*RecNum].timestamp_gps, &timestamp_gps_beg, sizeof(struct timeval));
      memcpy(&Rec[*RecNum].timestamp_pc, &timestamp_pc_beg, sizeof(struct timeval));
      Rec[*RecNum].AcqSeqNo = iteration_beg;
      Rec[*RecNum].blk_nano = blk_nano_beg;

      // fprintf(stdout,"===== BEAM >>>> RECNUM = %d BEG_SEQ = %d \n",*RecNum, Rec[*RecNum].AcqSeqNo);
      Rec[*RecNum].pc_time = timestamp_pc.tv_sec + timestamp_pc.tv_usec / 1e6;
      //  if(myrank == NUM_PROCESSES+2)       fprintf(pr_time, "PROCESS TIME = %lf\n", Rec[*RecNum].pc_time);
      set_rec_time(*RecNum, timestamp_pc);
      process(*RecNum);

      // After filling current Rec (RecNum) process that, unmark the previous Rec
      // -- YG/ABA 26/6/2001
      Rec[*RecNum - 1].rec_flag &= ~Marked;
      if (LogEvents)
        seq_log(*RecNum - 1);

      // if (*RecNum == (CurrRecs-1))
      if (*RecNum == (CurrRecs))
      {
        // fprintf(stderr, "## At RecNum= %d \n",  CurrRecs);
        // For RecNum = MaxRec, mark the 0th Rec ( Not available for read) then copy
        //       Rec of MaxRec(16th) to 0th. Rec of MaxRec is always marked.
        //       -- YG/ABA 26/6/2001.

        Rec[0].rec_flag = Rec[CurrRecs].rec_flag;
        //      set_blk_eqn();

        RackInfo.blocks = 0;
        RackInfo.time_err = 0;
        memmove((char *)Rec[0].begp - ExtraBuf,
                (char *)Rec[CurrRecs].begp - ExtraBuf, CurrRecSize + ExtraBuf);
        Rec[0].rec_time = Rec[CurrRecs].rec_time;
        Rec[0].rec_seq = Rec[CurrRecs].rec_seq;
        Rec[0].pc_time = Rec[CurrRecs].pc_time;

        memcpy(&Rec[0].timestamp_gps, &Rec[CurrRecs].timestamp_gps, sizeof(struct timeval));
        memcpy(&Rec[0].timestamp_pc, &Rec[CurrRecs].timestamp_pc, sizeof(struct timeval));
        Rec[0].AcqSeqNo = Rec[CurrRecs].AcqSeqNo;
        Rec[0].blk_nano = Rec[CurrRecs].blk_nano;

        *RecNum = 1;
        //     *RecNum = 0;
      }
      else
        ++*RecNum;
      Rec[*RecNum].rec_flag = Marked;
      Info->acq_flag &= ~Marked;
    }
    stoke_integ = 0;
    // if(strcmp(beam_mode, "IA") == 0) IAstoke_integ=0;
    // if(strcmp(beam_mode, "PA") == 0) PAstoke_integ=0;
    read_time = 1;
  }
  return 0;
}

// int initialize_IAPA(int Pols, int curr_chans, int curr_integ, int output_size, int ddc, int decimation, int beamid, nBMs)
int initialize_IAPA(int Pols, int curr_chans, int curr_integ, int output_size, int beamid, int nBMs)
{
  //int i, l, id_TEL_hdr = 0;
  int i, l;

  // id_TEL_hdr = shmget(BeamHdrKey, sizeof(BeamHeaderType), IPC_CREAT | 0666);
  // if (id_TEL_hdr < 0)
  // {
  //   perror("TEL shmget BEAM HEADER ID");
  //   return -1;
  // }
/*
  BeamHeaderType *dataHdr_TEL = (BeamHeaderType *)shmat(id_TEL_hdr, 0, 0);
  if (dataHdr_TEL == (void *)-1)
  {
    perror("TEL shmat BEAM HEADER");
    return -1;
  }
  dataHdr_TEL->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode = nBMs;
*/

  int Marker_Offsets[3][5] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
//  FILE *shmfile;
  long int CurrTotalWords = TotalWords;
  long int CurrShmDataSize = ShmDataSize;
  long int CurrShmSize = ShmSize;
//  int ref_beam_int = 128;
//  int ref_gpu_chans = 2048;

  // CurrRecs = Total No. of shm BLOCKS to be available in SHM.
  CurrRecs = 8; // 8 buffers of shared memory, hardcoded in acqpsr.h also.

  MaxTimeErr = CurrRecs / 2; // MaxTimeErr = 16/2
  // CurrSamples = (32*16*MAX_BLOCKS);  // MAX_BLOCKS = 8 no. of SHM BLOCKS, 32 factor in multiplication to match the beam data size, ad-hoc
  CurrSamples = (32 * 25); // 800 (no. of FFT in 1Sec.)

  // CurrTotalWords = (MaxPols * CurrSamples * ref_gpu_chans * Pols * ref_beam_int)/curr_integ; // independent of current channels.
  CurrTotalWords = CurrSamples * curr_chans; // independent of current channels.
  CurrRecSize = CurrTotalWords * WordSize / (output_size);
  CurrRecSize *= nBMs;
  fprintf(stdout, "## ----------------------------------------------- \n");
  // fprintf(stdout, "BEAM >>>> POLS=%d CURRSAMP = %d CHAN = %d CURRRECSIZE = %d CurrTotWords = %d \n", Pols, CurrSamples, curr_chans, CurrRecSize, CurrTotalWords);
  // fprintf(stdout, "BEAM >>>> OUT SIZE      = %d \n", output_size);
  // fprintf(stdout, "BEAM >>>> Ref CHANS     = %d \n", ref_gpu_chans);
  // fprintf(stdout, "BEAM >>>> MAXPOLS       = %d \n", MaxPols);
  //  fprintf(stdout, "BEAM >>>> NBEAMS        = %d \n", nBMs);
  // fprintf(stdout, "BEAM >>>> REF INTEG     = %d \n", ref_beam_int);
  fprintf(stderr, "# INTEG = %d, CHANS = %d DATA SIZE = %ld\n", curr_integ, curr_chans, CurrRecSize);
  //  fprintf(stdout, "BEAM >>>> INTEG         = %d \n", curr_integ);
  // fprintf(stdout, "BEAM >>>> POLS          = %d \n", Pols);
  // fprintf(stdout, "BEAM >>>> CURRSAMP      = %d \n", CurrSamples);
  // fprintf(stdout, "BEAM >>>> CHAN          = %d \n", curr_chans);
  // fprintf(stdout, "BEAM >>>> CurrTotWords  = %d \n", CurrTotalWords);
  //  fprintf(stdout, "BEAM >>>> CURRRECSIZE   = %d \n", CurrRecSize);
  fprintf(stdout, "## ----------------------------------------------- \n");

  CurrShmDataSize = (((CurrRecs + 1) * CurrRecSize + ExtraBuf) / PageSize + 1) * PageSize;
  CurrShmSize = ShmDataOff + CurrShmDataSize;

  // fprintf(stdout,"\nBEAM >>>> Shm information (acqpsr shm):\n");
  // fprintf(stdout,"BEAM >>>>   Shm Size     : %ld Bytes \n",CurrShmSize);
  // fprintf(stdout,"BEAM >>>>   Block Size   : %ld Bytes \n", CurrRecSize);
  // fprintf(stdout,"BEAM >>>>   Num of blocks: %d \n", CurrRecs);
  // fprintf(stdout,"BEAM >>>>   Time Samples : %d, curr_chans: %d, MaxPols: %d\n", CurrSamples, curr_chans, MaxPols);
  // system("ipcrm -M 0x0000040a");
  ShmId = shmget((int)(ShmKey + beamid), CurrShmSize, 0777 | IPC_CREAT);
  if (ShmId < 0)
  {
    perror("shmget failed:");
    return -1;
  }
  // write ShmId in the file
  //
  //  Not writing shm key in file, shm is not generatew with IPC_PRIVATE
  /*
    if((shmfile = file_fopen(IPCshare_file,"w","INITPATH")) == NULL)
    {
       perror("ShmIA file");
       exit(EXIT_FAILURE);
    }
    fprintf(shmfile,"ACQPSR_SHM %d\n",ShmId);
    fclose(shmfile);
  */

  Shmp = (char *)shmat(ShmId, 0, 0);
  if (Shmp == (void *)-1)
  {
    perror("shmat failed:");
    return -2;
  }
  if ((mlock(Shmp, CurrShmSize)) == -1)
  {
    perror("mlock failed\n");
  }

  /*
    if(errlog['p']==1)
    {
      fprintf(AcqErrLog,"  Shmp pointer =%p, ShmId=%d\n", Shmp, ShmId);
      fprintf(AcqErrLog,"  Shm key stored in file: %s \n", IPCshare_file);
    }
    fprintf(stdout,"  Shm pointer  : %p, Shm Id :%d\n", Shmp, ShmId);
    fprintf(stdout,"  Shm key stored in file    : %s \n\n", IPCshare_file );
  */

  // getchar();
  Info = (GlobalInfoType *)Shmp;
  Info->shmp = Shmp;
  Info->acq_flag = Marked | UnInitialized;
  Info->mark_num = 1;
  ISTeqn = Info->pc2ist;
  BLKeqn = Info->pc2blk;
  GPSval = Info->gps_val;
  GPSseq = Info->gps_seq;
  BLKval = Info->blk_val;
  BLKseq = Info->blk_seq;

  gettimeofday((struct timeval *)BlockInfo.time_val, (struct timezone *)TZ);
  TZoff = -TZ[0] * GPSmin;
  //  fprintf(stdout, "============> TZOFF = %ld\n", TZoff);
  //  fprintf(stdout, "============> TVAL0 = %ld\n", BlockInfo.time_val[0]);
  RefDayTime = BlockInfo.time_val[0] + TZoff -
               (BlockInfo.time_val[0] + TZoff) % 86400 - TZoff;
  //  fprintf(stdout, "============> RefDayTime = %ld\n", RefDayTime);
  /* TZoff is subtracted to include the TZ correction in RefDayTime */
  BlockInfo.time_val[0] -= RefDayTime;
  //  fprintf(stdout, "============> TVAL0 = %ld\n", BlockInfo.time_val[0]);
  Info->ref_time = RefDayTime;
  Info->tzoff = TZoff;

  if (init_gps() < 0)
    return -3; // init_gps does not do anything inside for moment, just return 0; // no GPS at moment.

  for (i = 0; i < MaxBLK; i++)
    BLKint[i].count = BLKseq[i] = 0;
  BLKtime = (double)FFT_CYCLE * DECIMATE * (CurrRecSize / (Pols * CHANNEL * 2)) / ((double)BASE_CLK); /*Nominal FFT_CYCLE=1024 BASE_CLK=33.33e6*/
  // BLKtime = (double)FFT_CYCLE*DECIMATE*Integ*(CurrRecSize/(Pols*CHANNEL*2))/((double)BASE_CLK); /*Nominal FFT_CYCLE=1024 BASE_CLK=33.33e6*/
  //  BLKtimeEps=(BLKtime*Pols)/(MaxPols*1.2);  // this will check if buffer is delayed by 0.209 sec (251mSec/1.2)
  BLKtimeEps = (2.5 * BLKtime * Pols) / (MaxPols); // this will check if buffer is delayed more than 2.5*251msec .
  BlktTimeEps = (2.5 * BLKtime * Pols) / (MaxPols);
  //  if(errlog['p']==1) fprintf(AcqErrLog,"BLKtime=%f <-some values are hardcoded in calc.\n", BLKtime);
  //  fprintf(stdout,"BLKtime=%f <-some values are hardcoded in calc.\n", BLKtime);
  // fprintf(stdout,"BEAM >>>> BLKtime=%f BLKtimeEps = %f.\n", BLKtime, BLKtimeEps);

  if (ISTeqnSet)
    Info->blk_time = BLKtime * ISTeqn[PC2IST];
  else
    Info->blk_time = BLKtime;

  Info->rec_ind = 0;
  RecNum = &Info->rec_ind;
  Rec = Info->rec;
  for (i = 0; i <= CurrRecs; i++)
  {
    Rec[i].rec_flag = UnInitialized | Marked;
    Rec[i].rec_seq = 0;
    Rec[i].beg_off = ShmDataOff + i * CurrRecSize;
    Rec[i].begp = (unsigned short *)(Shmp + Rec[i].beg_off);
    Rec[i].pc_time = 0;
    Rec[i].rec_time = 0;
  }

  RackInfo.data_flag = 0;
  RackInfo.dbuf_seq = 0;
  RackInfo.blocks = 0;
  RackInfo.time_err = 0;
  RackInfo.prev_time = 0;

  BlockTime = RackInfo.block_time.tv;
  BlockTime[0] = BlockInfo.time_val[0];
  BlockTime[1] = BlockInfo.time_val[1];

  for (i = 0; i < 3; i++)
    for (l = 0; l < 5; l++)
      Info->marker_offsets[i][l] = Marker_Offsets[i][l];

  Info->acq_flag &= ~UnInitialized;

  return ShmId;
}

/*****************************************************************
 * FUNCTION: void set_blk_times(double pc_time, int seq, int time_ok)
 *
 ******************************************************************/
int set_seq_num(int rec_num, double pc_time)
{
  int ind, seq, time_ok = 0;
  double diff, rseq;
  static unsigned long rseq_new = 0;

  ind = BLKind ? BLKind - 1 : MaxBLK - 1;
  diff = fabs(pc_time - BLKval[ind]);

  rseq = ToBLK(pc_time);
  seq = (int)(rseq + 0.5);
  //  Rec[rec_num].rec_seq = rseq_new++;
  Rec[rec_num].rec_seq = rseq_new;
  // fprintf(stdout, "BEAM >>>> RECORD SEQUENCE = %ld Block_ind = %d\n", rseq_new, Info->blk_ind); // SSK COMMENTED
  rseq_new++;
  return 0;

  fprintf(stdout, "BEAM >>>> ======>>>>> IND = %2d BLKSEQ[IND] = %d RSEQ=%ld\n", ind, BLKseq[ind], rseq_new);

  if (fabs(rseq - seq) > BLKtimeEps / BLKtime)
  {
    fprintf(stdout, "BEAM >>>> Large diff. in BLKseq num, %d %f, rec=%d pct=%f\n", seq, rseq, rec_num, pc_time);

    //    if (RackInfo.time_err >= MaxTimeErr)
    //    { fprintf(FE,"prev_time changed from %f to %f\n",
    //              RackInfo.prev_time, ToPCT(seq+1));
    //      RackInfo.prev_time = ToPCT(seq+1);
    //    }
  }
  else
    time_ok = 1;
  time_ok = 1; // All times are OK
  if (abs(seq - BLKseq[ind]) > MaxBLK / 2)
  {
    fprintf(stdout, "BEAM >>>> A big break in data seq, seq=%d, prev_seq=%d, pc_time=%f\n",
            seq, BLKseq[ind], pc_time);
  }
  if (diff > BLKtime / 2)
    set_blk_times(pc_time, seq, time_ok); // original
  return 0;
}

/*****************************************************************
 * FUNCTION: int set_rec_time(int rec_num, struct timeval timestamp)
 *
 ******************************************************************/
int set_rec_time(int rec_num, struct timeval timestamp)
{
  int blocks;
  double curr_time, exp_time;

  get_block_info(&BlockInfo);
  RackInfo.data_flag = 0;
  blocks = BlockInfo.count - RackInfo.dbuf_seq;
  //  if(myrank == NUM_PROCESSES+2) fprintf(stdout, "================= >>>> BLOKCS = %d BLOCKSPERREC = %d\n",blocks,BlocksPerRec);
  RackInfo.dbuf_seq = BlockInfo.count;
  //  memcpy(BlockTime, BlockInfo.time_val, 2*sizeof(long));
  memcpy(BlockTime, &timestamp, sizeof(struct timeval));
  BlockTime[0] -= RefDayTime;
  curr_time = BlockTime[0] + BlockTime[1] / 1e6;

  /* Previous Buffer Time is supplied as the current buffer time because
     the device driver supplies the time that corresponds to the end of
     arrival of the data blocks at the card. In our case, for the last
     block of the BLK cycle, it corresponds to the end of BLK cycle.
  */
  Rec[rec_num].pc_time = RackInfo.prev_time;
  exp_time = RackInfo.prev_time + (BLKtime * blocks) / BlocksPerRec;
  //  fprintf (stdout, "=====>>>>>>>> BLKtime = %f BLOCKS = %d  BLOCKSPERRCS = %d\n", BLKtime, blocks, BlocksPerRec);
  //  if(myrank == NUM_PROCESSES+2)  fprintf (stdout, "=====SSSSSSSS MYRANK = %d BLOCKS = %d  BLOCKSPERRCS = %d\n", myrank, blocks, BlocksPerRec);
  if (blocks != BlocksPerRec)
  {
    RackInfo.data_flag |= BlockErr;
    //
    //   if(errlog['p']==1) {
    //     fprintf(AcqErrLog,"Wrong Blocks=%d, Expected=%d\n", blocks, BlocksPerRec);
    //     fprintf(AcqErrLog,"  dbuf_seq=%d, curr=%12.6f, exp=%12.6f\n", RackInfo.dbuf_seq, curr_time, exp_time);
    //   }
    //
    fprintf(stdout, "BEAM >>>> Wrong Blocks=%d, Expected=%d\n", blocks, BlocksPerRec);
    fprintf(stdout, "BEAM >>>>   dbuf_seq=%d, curr=%12.6f, exp=%12.6f diff=%f\n", RackInfo.dbuf_seq, curr_time, exp_time, fabs(curr_time - exp_time));
    // Considering that the wrong blocks may have been reported due
    // to the spurious interrupt from HSDIO. Accept curr_time as right.
    // exp_time = curr_time;
    // RackInfo.time_err++;
    //
  }
  if (fabs(curr_time - exp_time) > BLKtimeEps)
  {
    //    if(errlog['p']==1) fprintf(AcqErrLog,"BlkTErr dbuf_seq=%d, curr=%12.6f, exp=%12.6f\n", RackInfo.dbuf_seq, curr_time, exp_time);
    //  if(myrank == NUM_PROCESSES+2) fprintf(stdout,"BLKtimeEps = %f BlkTErr dbuf_seq=%d, curr=%12.6f, exp=%12.6f DIFF = %f\n", BLKtimeEps, RackInfo.dbuf_seq, curr_time, exp_time, curr_time-exp_time);
    RackInfo.data_flag |= TimeErr;
    if (++RackInfo.time_err >= MaxTimeErr)
      RackInfo.prev_time = curr_time;
    else
      RackInfo.prev_time = exp_time;
  }
  else
    RackInfo.prev_time = curr_time;
  //  if(myrank == NUM_PROCESSES+2)	fprintf(stdout, "==========>>>> SET_REC_TIME ===  RACKINFO.BLOKCS  = %d BLOCKS = %d\n", RackInfo.blocks, blocks);
  RackInfo.blocks += blocks;
  return 0;
}

/*****************************************************************
 * FUNCTION: int get_rec(int rec_num, int accept, struct timeval timestamp)
 *
 ******************************************************************/
int get_rec(int rec_num, int accept, struct timeval timestamp)
{
  //int n = 0, m = 0;
  // data is already read, at a place where this routine is called from.
  if (accept)
  {
    double curr_time;
    get_block_info(&BlockInfo);
    RackInfo.dbuf_seq = BlockInfo.count;
    curr_time = BlockInfo.time_val[0] - RefDayTime +
                BlockInfo.time_val[1] / 1e6;

    //  start time of data block
    Rec[rec_num].pc_time = curr_time - BLKtime;
    RackInfo.prev_time = curr_time;
    return 0;
  }
  return set_rec_time(rec_num, timestamp);
}

/*****************************************************************
 * FUNCTION: int process(int rec_num)
 *
 ******************************************************************/

int process(int rec_num)
{
  RecType *rec = &Rec[rec_num];

  /*
     Time stamped here is a true IST when the first word of the
     first BLK cycle landed at FFT. This time is close to few 10s
     of microsecond after the signal arrived at an antenna =>
     signal carrier delay + Delay in DPC + delay of processing
     between DPC and Pulsar Receiver.
  */
  rec->rec_flag |= RackInfo.data_flag | GoodData;
  rec->rec_time = ToIST(rec->pc_time);
  set_seq_num(rec_num, rec->pc_time);
  return 0;
}

/*****************************************************************
 * FUNCTION: void set_blk_times(double pc_time, int seq, int time_ok)
 *
 ******************************************************************/
void set_blk_times(double pc_time, int seq, int time_ok)
{
  int i, l, n, ind;
  ind = BLKind ? BLKind - 1 : MaxBLK - 1;
  if ((n = seq - BLKseq[ind]) <= 0)
    return;
  l = BLKseq[ind];
  ind = BLKind;
  for (i = 0; i < n; i++)
  {
    BLKseq[ind] = ++l;
    if (i == n - 1 && time_ok)
      BLKval[ind] = pc_time;
    else
      BLKval[ind] = ToPCT(l);
    if (++ind >= MaxBLK)
      ind = 0;
  }
  Info->blk_ind = BLKind = ind;
}

/*****************************************************************
 * FUNCTION: void init_blk_eqn(int rec_num)
 *
 ******************************************************************/
void init_blk_eqn(int rec_num)
{
  BLKeqn[PC2BLK] = 1 / BLKtime;
  BLKeqn[BLKref] = DataSeq;
  BLKeqn[PCref] = Rec[rec_num].pc_time;
  // BLKeqnSet      = 1;

  BLKind = 0;
  BLKval[BLKind] = BLKeqn[PCref];
  BLKseq[BLKind] = BLKeqn[BLKref];

  Info->blk_ind = BLKind++;
  Info->blk_time = ISTeqn[PC2IST] * BLKtime;
  //  fprintf(stdout,"BEAM >>>> BLK time: PC Clock=%9.6f, IST=%9.6f\n", BLKtime, Info->blk_time);
}

/*****************************************************************
 * FUNCTION: void get_intr_times(void)
 *
 ******************************************************************/
void get_intr_times(void)
{
  double curr_time;
  static double max_time = 0, last_time = 0;
  struct timeval tv;
//  int intr = 0;

  GetIntr = 0;
  gettimeofday(&tv, NULL);
  curr_time = tv.tv_sec - RefDayTime + tv.tv_usec / 1e6;
  /*
    while((curr_time - last_time) < BLKtime) {
                  usleep(500);
            gettimeofday(&tv, NULL);
            curr_time = tv.tv_sec - RefDayTime + tv.tv_usec/1e6;
          }
  */
  if (last_time > 0 && curr_time - last_time > max_time)
  {
    max_time = curr_time - last_time;
    //    if(errlog['p']==1) fprintf(AcqErrLog,"get_intr is called after %10.6f sec\n", max_time);
    fprintf(stdout, "BEAM >>>> get_intr is called after %10.6f sec\n", max_time);
  }
  last_time = curr_time;
  /*
    if (GPS)
    { static double last_gps=0;
      if (last_gps == 0) // first time
        last_gps = GPSval[GPSind ? GPSind-1 : MaxGPS-1];
      if (curr_time - last_gps - 2*TickTime > GPStime)
      { if ((intr = get_gps_times(curr_time)) < 0) NotStop = 0;
        else if (intr > 0)
        { last_gps = GPSval[GPSind ? GPSind-1 : MaxGPS-1]; set_ist_eqn(); }
      }
    }
  */
}

/*****************************************************************
 * FUNCTION: int set_blk_eqn(void)
 *
 ******************************************************************/
int set_blk_eqn(void)
{
  static int blk_err = 0;
  double period, ref_val, rseq;
  int m, n, seq;

  if (RackInfo.blocks != BlocksPerRec * CurrRecs || RackInfo.time_err >= MaxTimeErr)
  {
    //    if(errlog['p']==1) fprintf(AcqErrLog,"BLKeqn not set, dbuf_seq=%d time_errs=%d" " block_errs=%d\n", RackInfo.dbuf_seq, RackInfo.time_err, RackInfo.blocks);
    fprintf(stdout, "BEAM >>>> BLKeqn not set, dbuf_seq=%d time_errs=%d"
                    " block_errs=%d\n",
            RackInfo.dbuf_seq,
            RackInfo.time_err, RackInfo.blocks);
    return -1;
  }
  else
  {
    for (m = 0; m <= CurrRecs; m++)
      if ((Rec[m].rec_flag & (TimeErr | BlockErr)) == 0)
        break;
    for (n = CurrRecs; n > m; n--)
      if ((Rec[n].rec_flag & (TimeErr | BlockErr)) == 0)
        break;
    if (n > m)
      period = (Rec[n].pc_time - Rec[m].pc_time) / (n - m);
    else
      return -2;
  }
  if (fabs(period - BLKtime) > BlktTimeEps)
  {
    blk_err++;
    //    if(errlog['p']==1) fprintf(AcqErrLog,"BLKeqn not set, dbuf_seq=%d blk_err=%d period=%9.6f" " BLKtime=%9.6f\n", RackInfo.dbuf_seq, blk_err, period, BLKtime);
    fprintf(stdout, "BEAM >>>> BLKeqn not set, dbuf_seq=%d blk_err=%d period=%9.6f"
                    " BLKtime=%9.6f\n",
            RackInfo.dbuf_seq,
            blk_err, period, BLKtime);
    {
      enum
      {
        Max = 10
      };
      int i, j, k;
      static int ind = 0;
      static double blk_time[Max] = {0};
      blk_time[ind] = period;
      if (++ind == Max)
        ind = 0;
      if (blk_err >= Max)
      {
        double p = 0.0;
        for (i = ind; i < ind + Max / 2; i++)
        {
          k = 0;
          p = blk_time[i % Max];
          for (j = i + 1; j < ind + Max; j++)
            if (fabs(blk_time[i % Max] - blk_time[j % Max]) < BlktTimeEps)
            {
              k++;
              p += blk_time[j % Max];
            }
          if (k >= Max / 2)
          {
            p /= k;
            break;
          }
        }
        if (i >= ind + Max / 2)
          return -2;
        else
          period = p;
      }
    }
  }

  ref_val = Rec[n].pc_time;

  rseq = ToBLK(ref_val);
  seq = (int)(rseq + 0.5);
  if (rseq - seq > BLKtimeEps / BLKtime)
  {
    //    if(errlog['p']==1) fprintf(AcqErrLog,"Large change in blk_eqn_seq (%f, %d)\n", rseq, seq);
    fprintf(stdout, "BEAM >>>> Large change in blk_eqn_seq (%f, %d)\n", rseq, seq);
  }
  BLKtime = (99 * BLKtime + period) / 100;
  Info->blk_time = ISTeqn[PC2IST] * BLKtime;
  BLKeqn[PC2BLK] = 1 / BLKtime;
  BLKeqn[PCref] = ref_val;
  BLKeqn[BLKref] = seq;
  if (LogEvents)
    blk_log();
  /* fprintf(FE,"Next BLK eqn: PC2BLK=%f, PCref=%f, BLKref=%f\n",
              BLKeqn[PC2BLK], BLKeqn[PCref], BLKeqn[BLKref]); $$ */
  blk_err = 0;
  return 0;
}

/*****************************************************************
 * FUNCTION:
 *
 *  write BLKlog to file LogBuf
 ******************************************************************/
void blk_log(void)
{
  *(double *)BLKlog = BLKeqn[PCref];
  BLKlog->flag = BLKTag;
  BLKlog->blk_ref = (unsigned)BLKeqn[BLKref];
  BLKlog->pc2blk = (float)BLKeqn[PC2BLK];
  //if (fwrite(LogBuf, sizeof(LogBuf), 1, Logfd) < 0)
  //{ perror("fwrite BLK LogBuf"); LogEvents = 0; }
}

/*****************************************************************
 * FUNCTION: void seq_log(int rec_num)
 *
 *  write SeqLog to file LogBuf
 ******************************************************************/
void seq_log(int rec_num)
{
  SeqLog->dummy = 0;
  SeqLog->rec_num = rec_num;
  *(double *)SeqLog = Rec[rec_num].rec_time;
  SeqLog->flag = (Rec[rec_num].rec_flag & ~FlagMask) | SeqTag;
  SeqLog->seq_num = Rec[rec_num].rec_seq;
  //  if (fwrite(LogBuf, sizeof(LogBuf), 1, Logfd) < 0) { perror("fwrite Seq LogBuf"); LogEvents = 0; }
}

int init_gps(void)
{
  int i, l, ref = 0;

  //  return 0; // later any init_gps routines can be coded. This routine is from existing
  // hardware correlator pulsar pipeline. GSB may not be keeping the track of GPS
  // the same way as GHB. Need to think.
  //
  for (i = 0; i < MaxGPS; i++)
    GPSint[i].count = GPSseq[i] = 0;
  if (GPS)
  {
    GPSfd = open(GPSdev, O_RDONLY);
    fprintf(stdout, "BEAM >>>> ==============> INSIDE INIT_GPS\n");
    if (GPSfd < 0)
      perror(GPSdev);
    GPS = GPSfd >= 0;
  }
  if (GPS)
  { /* Do not proceed without having at least three GPS min pulses */
#ifdef _SIM_GPS_
    {
      int arg[3] = {5, (long)(GPSmin * 1e6), (long)(GPStimeEps * 1e6)};
      /* 5 pulses before */
      ioctl(GPSfd, SetTime + (3 * sizeof(int) << 16), arg);
    }
#endif
    Info->acq_flag |= GPSpresent;
    Gind = 0;
    GPSind = 0;
    GPStime = 0;
    l = read(GPSfd, GPSint, (MaxGPS - 2) * TimeSz);
    if (l > 0)
      Gind = l / TimeSz;
    for (i = 0; i < Gind; i++)
      GPSival[i] = GPSval[i] =
          (GPSint[i].tv[0] - RefDayTime) + GPSint[i].tv[1] / 1e6;
    do
    {
      double gps_time = 0;
      if (Gind >= 3)
      {
        for (i = 2; i < Gind; i++)
          if ((gps_time = test_gps(i)) > 0)
            break;
        if (i < Gind)
        {
          GPStime = gps_time;
          ref = i - 2;
          break;
        }
      }
      while (1)
      {
        l = read(GPSfd, GPSint + Gind, TimeSz);
        if (l <= 0 || (Gind > 0 && GPSint[Gind].seq == GPSint[Gind - 1].seq))
          usleep(1000000);
        else
          break;
      }
      GPSival[Gind] = GPSval[Gind] =
          (GPSint[Gind].tv[0] -= RefDayTime) + GPSint[Gind].tv[1] / 1e6;
      Gind++;
    } while (Gind < MaxGPS - 2);

    fprintf(stdout, "BEAM >>>> Initial GPS times are: ref=%d\n", ref);
    for (i = 0; i < Gind; i++)
    {
      fprintf(stdout, "BEAM >>>> %3d %7d %12.6f\n", i, GPSint[i].seq, GPSival[i]);
    }
    if (GPStime == 0)
    {
      fprintf(stdout, "BEAM >>>> GPS times are not blkble within %f sec\n", BlktTimeEps);
      return -1;
    }
    else if (fabs(GPSmin - GPStime) > GPStimeEps)
    {
      fprintf(stdout, "BEAM >>>> GPStime=%10.6f sec; Bad PC Clock\n", GPStime);
      return -2;
    }
    else
    {
      fprintf(stdout, "BEAM >>>> GPStime determined to be %10.6f\n", GPStime);
    }

    if (ref > 0)
    {
      for (i = ref; i < Gind; i++)
      {
        GPSint[i - ref] = GPSint[i];
        GPSival[i - ref] = GPSval[i - ref] = GPSival[i];
      }
      Gind -= ref;
    }
    GPSind = Gind;
    /* GPSval and GPSival are treated differently from here */
    for (i = 3; i < GPSind;)
    {
      l = set_gps(i);
      if (l < 0 || GPSind == 0)
      {
        fprintf(stdout, "BEAM >>>> Multiple GPS timer error during initialization\n");
        return -3;
      }
      i = (i + l) % MaxGPS;
    }
    /* GPSseq is deliberately set after the set_gps, set_gps uses GPSseq */
    GPSseq[GPSind - 1] = (int)((GPSval[GPSind - 1] + GPSmin / 2) / GPSmin);
    for (i = GPSind - 2; i >= 0; i--)
      GPSseq[i] = GPSseq[i + 1] - 1;

    ISTeqn[PC2IST] = GPSmin / GPStime;
    ISTeqn[PCref] = GPSval[GPSind - 1];
    ISTeqn[ISTref] = GPSseq[GPSind - 1] * GPSmin;

    if (GPSind == MaxGPS)
      GPSind = 0;
    Info->gps_ind = GPSind;
    if (LogEvents)
      ist_log();
  }
  else
  {
    ISTeqn[PC2IST] = 1;
    ISTeqn[PCref] = 0;
    ISTeqn[ISTref] = 0;
  }
  ISTeqnSet = 1;
  // fprintf(stdout,"BEAM >>>> Initial IST eqn: PC2IST=%f, PCref=%f, ISTref=%f\n", ISTeqn[PC2IST], ISTeqn[PCref], ISTeqn[ISTref]);
  return 0;
}

double test_gps(int ind2)
{
  int ind0, ind1;
  ind0 = (ind2 - 2 + MaxGPS) % MaxGPS;
  ind1 = (ind2 - 1 + MaxGPS) % MaxGPS;
  if (fabs(2 * GPSival[ind1] - (GPSival[ind0] + GPSival[ind2])) < BlktTimeEps &&
      fabs(GPSival[ind2] - GPSival[ind1] - GPSmin) < GPStimeEps)
    return (GPSival[ind2] - GPSival[ind0]) / 2;
  return 0;
}

int set_gps(int i)
{
  int j;
  double diff, val;
  if ((j = i - 1) < 0)
    j = MaxGPS - 1;
  val = GPSval[j] + GPStime;
  if (fabs(GPSval[i] - val) > GPStimeEps)
  {
    diff = GPSval[i] - GPSval[j];
    fprintf(stdout, "BEAM >>>> Bad GPS time, GPSval[%2d]=%10.6f, GPSval[%2d]=%10.6f,"
                    " diff=%10.6f, GPSseq[%2d]=%d\n",
            j, GPSval[j],
            i, GPSval[i], diff, i, GPSseq[i]);
    if (diff < GPStime - 1.5 * TickTime)
    {
      fprintf(stdout, "BEAM >>>>   Discarded GPSval=%10.6f, GPSseq=%d\n", GPSval[i], GPSseq[i]);
      shift(i, 1, GPSind % MaxGPS, MaxGPS, GPSval);
      if (GPSind)
        --GPSind;
      else
        GPSind = MaxGPS - 1;
      GPSseq[GPSind] = 0;
      return 0;
    }
    else
    {
      if (ISTeqnSet)
        val = ToPC(GPSseq[i] * 60);
      if (diff > GPStime + 1.5 * TickTime)
      {
        fprintf(stdout, "BEAM >>>>   Added GPSval=%10.6f, GPSseq=%d\n", val, GPSseq[i]);
        shift(i, -1, GPSind % MaxGPS, MaxGPS, GPSval);
        GPSseq[BLKind] = GPSseq[i] + 1;
        if (++GPSind == MaxGPS)
          GPSind = 0;
      }
      else
      {
        fprintf(stdout, "BEAM >>>>   GPSval[%d]=%12.6f => %12.6f GPSseq[%d]=%d\n",
                i, GPSval[i], val, i, GPSseq[i]);
      }
      GPSval[i] = val;
    }
  }
  return 1;
}

void shift(int ind, int step, int last, int max, double *val)
{
  int from, to;
  if (step > 0)
  {
    to = ind, from = (ind + step) % max;
    while (from != last)
    {
      val[to] = val[from];
      if (++from == max)
        from = 0;
      if (++to == max)
        to = 0;
    }
  }
  else if (step < 0)
  {
    if ((from = last - 1) < 0)
      from = max - 1;
    to = (from - step) % max;
    do
    {
      val[to] = val[from];
      if (from == ind)
        break;
      if (--from < 0)
        from = max - 1;
      if (--to < 0)
        to = max - 1;
    } while (1);
  }
}

void ist_log(void)
{
  *(double *)ISTlog = ISTeqn[PCref];
  ISTlog->flag = ISTTag;
  ISTlog->ist_ref = (unsigned)ISTeqn[ISTref];
  ISTlog->pc2ist = (float)ISTeqn[PC2IST];
  //if (fwrite(LogBuf, sizeof(LogBuf), 1, Logfd) < 0)
  //{ perror("fwrite IST LogBuf"); LogEvents = 0; }
}
//========================= END of ACQPSR ROUTINES  ================================//
