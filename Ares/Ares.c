/******************************************************************************************************************************************************
 * Ares.c
 *
 * Software to read beam data from SPOTLIGHT multi-beam correlator (TEL_SHM) and create the FRB_SHM.
 * COMPILE: gcc -g -O3 -Wall -Wextra -std=gnu99 Ares.c -o Ares -lm -lc
 * USAGE: ./Ares [observation duration (seconds) = 3600]
 * Author: Kenil Ajudiya. Contact: kenilr@iisc.ac.in                                    Date: July 27, 2024
 ******************************************************************************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/shm.h>
#include "gmrt_newcorr.h"
#include "acqpsr.h"
#include "frb_shm.h"

#define BeamHdrKey 1050

static BeamHeaderType *dataHdr_TEL;
static RecType *Rec;
GlobalInfoType *dataBuffer_TEL;
BeamHeaderType *dataHdr_FRB;
DataBuffer *dataBuffer_FRB;
unsigned char *Shmp, *blkp_TEL, *blkp_FRB;
unsigned int RecNum = 0, TEL_to_FRB_RecNum = 0, RecNum_FRB = 0, DataSeq = 0, DataSeq_FRB = 0;
// RecNum goes from 0 to MaxRecs (MaxRecs = 8 for TEL_SHM and 12*32 for FRB_SHM.)
// Only for the sake of this program, I am pretending as if FRB_SHM is 12 * 32 blocks long, with each block equal in size to that of TEL_SHM.
// In the FRB_SHM, the curRecord and curBlock are updated after accumulating 32 blocks of TEL_SHM.
// DataSeq keeps track of the number of blocks (or records) read since the starting of the observation. It is, in a sense, cummulative RecNum.

unsigned int *shm_marker_p;
double Data_Time[MaxRecs + 1] = {0.0};
int ind = 0, beam_ID = 0;
struct tm *local_t;
char time_string[40];
int usec = 0;
double blk_nano = 0;

struct timeval tv1;
struct timezone tz1;
int nBeams;

int initialize_TEL_SHM(void)
{
  int id_TEL_hdr = 0, id_TEL_buf = 0;

  id_TEL_hdr = shmget(BeamHdrKey, sizeof(BeamHeaderType), SHM_RDONLY);
  if (id_TEL_hdr < 0)
  {
    perror("TEL shmget BEAM HEADER ID");
    return -1;
  }
  dataHdr_TEL = (BeamHeaderType *)shmat(id_TEL_hdr, 0, 0);
  if (dataHdr_TEL == (void *)-1)
  {
    perror("TEL shmat BEAM HEADER");
    return -1;
  }

  nBeams = dataHdr_TEL->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode;

  int local_ExtraBuf = 64;
  long CurrShmDataSize = (((MaxRecs + 1) * TEL_SHM_BlockSize * nBeams + local_ExtraBuf) / PageSize + 1) * PageSize;
  int local_ShmDataOff = ((sizeof(GlobalInfoType) + local_ExtraBuf) / PageSize + 1) * PageSize;
  long CurrShmSize = local_ShmDataOff + CurrShmDataSize;

  id_TEL_buf = shmget(ShmKey, CurrShmSize, 0);
  if (id_TEL_hdr < 0)
  {
    perror("TEL shmget BEAM BUFFER ID");
    return -1;
  }
  dataBuffer_TEL = (GlobalInfoType *)shmat(id_TEL_buf, 0, SHM_RDONLY);
  if ((void *)dataBuffer_TEL == (void *)-1)
  {
    perror("TEL shmat BEAM BUFFER ID");
    return -3;
  }
  Shmp = (unsigned char *)dataBuffer_TEL;

  Rec = dataBuffer_TEL->rec;
  fprintf(stderr, "dataBuffer_TEL->rec_ind: %d\n", dataBuffer_TEL->rec_ind);
  RecNum = ((unsigned int)(dataBuffer_TEL->rec_ind - 1)) % MaxRecs;
  fprintf(stderr, "RecNum: %u\n", RecNum);
  DataSeq = Rec[RecNum].rec_seq; /* expected first rec */

  return 0;
}

int initialize_FRB_SHM()
{
  int id_FRB_hdr, id_FRB_buf;
  id_FRB_hdr = shmget(DasHeaderKey, sizeof(BeamHeaderType), IPC_CREAT | 0666);
  if (id_FRB_hdr < 0)
  {
    perror("FRB shmget HEADER ID");
    return -1;
  }
  dataHdr_FRB = (BeamHeaderType *)shmat(id_FRB_hdr, 0, 0);
  if ((void *)dataHdr_FRB == (void *)-1)
  {
    perror("FRB shmat HEADER ID");
    return -3;
  }

  id_FRB_buf = shmget(DasBufferKey, sizeof(DataBuffer), IPC_CREAT | 0666);
  if (id_FRB_hdr < 0)
  {
    perror("FRB shmget BUFFER ID");
    return -1;
  }
  dataBuffer_FRB = (DataBuffer *)shmat(id_FRB_buf, 0, 0);
  if ((void *)dataBuffer_FRB == (void *)-1)
  {
    perror("FRB shmat BUFFER ID");
    return -3;
  }

  // Copying the TEL_SHM header as it is for right now. This needs to be updated multiple times during the observation.
  memcpy(dataHdr_FRB, dataHdr_TEL, sizeof(BeamHeaderType));

  dataBuffer_FRB->timestamp_gps[0] = Rec[RecNum].timestamp_gps;
  dataBuffer_FRB->blk_nano[0] = Rec[RecNum].blk_nano; // Adding nanoseconds accuracy for the timestamps.
  dataBuffer_FRB->pc_time = Rec[RecNum].pc_time;
  dataBuffer_FRB->ref_time = dataBuffer_TEL->ref_time;
  dataBuffer_FRB->rec_time = Rec[RecNum].rec_time;

	dataBuffer_FRB->nBeams = nBeams;
  dataBuffer_FRB->curRecord = 0;
  dataBuffer_FRB->curBlock = 0;
  dataBuffer_FRB->blockSize = DataSize;
  dataBuffer_FRB->active = 1;
  dataBuffer_FRB->is_buf_empty = 1; // There is no meaningful data in the buffer. After filling in one block of buffer, this flag will be set to 0.
  dataBuffer_FRB->flag = 0;

  return 0;
}

int process_acq(void)
{
  if (dataBuffer_TEL->acq_flag & AcqOver) // Ask Sanjay whether this is useful here.
  {
    fprintf(stderr, "ACQOVER\n");
    return -2;
  }

  ind = dataBuffer_TEL->blk_ind > 0 ? dataBuffer_TEL->blk_ind - 1 : MaxBLK - 1;
  if (dataBuffer_TEL->blk_seq[ind] > DataSeq + MaxRecs - 1)
  {
    fprintf(stderr, "Processing lagged behind: DataSeq=%d AcqSeq=%d, RecNum=%d.\n", DataSeq, dataBuffer_TEL->blk_seq[ind], RecNum);
    fprintf(stderr, "New Seq. for processing=%d\a\n", dataBuffer_TEL->blk_seq[ind] - MaxRecs + 1);
    DataSeq = dataBuffer_TEL->blk_seq[ind] - MaxRecs + 1;
    RecNum = (dataBuffer_TEL->rec_ind - MaxRecs + 1) % MaxRecs;
  }

  while (1)
  {
    if (Rec[RecNum].rec_flag & Marked)
    {
      return -3;
    }
    if (Rec[RecNum].rec_seq < DataSeq)
    {
      RecNum = (RecNum + 1) % MaxRecs;
    }
    else
    {
      break;
    }
  }

  if (Rec[RecNum].rec_seq > DataSeq)
  {
    DataSeq++;
    return -4;
  }

  if ((Rec[RecNum].rec_flag & GoodData) == 0)
  {
    fprintf(stderr, "Bad data: DataSeq=%d flag=%4d\n", DataSeq, Rec[RecNum].rec_flag);
    ++DataSeq;
    return -5;
  }

  blkp_TEL = Shmp + Rec[RecNum].beg_off;
  Data_Time[RecNum] = Rec[RecNum].pc_time;
  local_t = localtime(&Rec[RecNum].timestamp_gps.tv_sec);
  gettimeofday(&tv1, &tz1);
  Data_Time[RecNum] = ((tv1.tv_sec - tz1.tz_minuteswest * 60) % 86400 + tv1.tv_usec / 1e6);

  if (TEL_to_FRB_RecNum == 0 && dataBuffer_FRB->flag != 0) // Data for 1 block of FRB_SHM has accumulated. Update the curRecord and curBlock.
  {
    DataSeq_FRB++;
    RecNum_FRB = (RecNum_FRB + 1) % MaxDataBlocks;

    dataBuffer_FRB->timestamp_gps[RecNum_FRB] = Rec[RecNum].timestamp_gps;
    dataBuffer_FRB->blk_nano[RecNum_FRB] = Rec[RecNum].blk_nano; // Adding nanoseconds accuracy for the timestamps.
    dataBuffer_FRB->pc_time = Rec[RecNum].pc_time;
    dataBuffer_FRB->ref_time = dataBuffer_TEL->ref_time;
    dataBuffer_FRB->rec_time = Rec[RecNum].rec_time;

    dataBuffer_FRB->curRecord = RecNum_FRB;
    dataBuffer_FRB->curBlock = DataSeq_FRB;
    dataBuffer_FRB->is_buf_empty = 0; // Now it is sure that there is at least one block of buffer filled with data.

    fprintf(stderr, "\n############ Block number in FRB_SHM: %d ############", dataBuffer_FRB->curRecord);
    fprintf(stderr, "\n############    Total blocks read: %d    ############\n", dataBuffer_FRB->curBlock);
  }
  // Each block of TEL_SHM will be copied as soon as it is available in the SHM.
  // The curRecord and other parameters are updated (see the above if condition) only after accumulating 32 blocks of TEL_SHM.
  // Pretending as if FRB_SHM is 12 * 32 blocks long, with each block equal in size to that of TEL_SHM.

  for (beam_ID = 0; beam_ID < nBeams; beam_ID++)
  {
    blkp_FRB = dataBuffer_FRB->data + ((long)DataSize * nBeams * RecNum_FRB) + ((long)DataSize * beam_ID) + ((long)TEL_SHM_BlockSize * TEL_to_FRB_RecNum);
    memcpy(blkp_FRB, blkp_TEL + (long)TEL_SHM_BlockSize * beam_ID, TEL_SHM_BlockSize);
  }

  usec = (Rec[RecNum].timestamp_gps.tv_usec);
  blk_nano = (Rec[RecNum].blk_nano);
  strftime(time_string, sizeof(time_string), "%Y %m %d %H %M %S", local_t);
  fprintf(stderr, "# LOG # Current data timestamp = %02d:%02d:%02d.%06d%03d \t Rec[RecNum].rec_seq: %u \t RecNum: %u\n", local_t->tm_hour, local_t->tm_min, local_t->tm_sec, usec, (int)(1 * blk_nano), Rec[RecNum].rec_seq, RecNum);

  DataSeq++;
  RecNum = (RecNum + 1) % MaxRecs;
  TEL_to_FRB_RecNum = (TEL_to_FRB_RecNum + 1) % TEL_to_FRB_Block_factor;
  dataBuffer_FRB->flag = 1;
  return 0;
}

int main(int argc, char *argv[])
{
  long obs_duration = round(3600 / 1.048576);
  if (argc != 2)
  {
    fprintf(stderr, "\nUSAGE: %s [observation duration (seconds) ]\n\n", argv[0]);
    fprintf(stderr, "Reading 1 hour of data by default.\n\n");
    return (-1);
  }
  else
    obs_duration = round(atoi(argv[1]) / 1.048576);

  if (initialize_TEL_SHM() < 0)
  {
    fprintf(stderr, "TEL_SHM initialization failure\n");
    return -3;
  }

  if (initialize_FRB_SHM() < 0)
  {
    fprintf(stderr, "FRB_SHM initialization failure\n");
    return -3;
  }

  fprintf(stderr, "-------------------------------------------------------------------------------------\n");
  fprintf(stderr, "   Num. Beams = %d TEL_SHM BLOCK SIZE = %ld, FRB_SHM BLOCK SIZE = %ld\n", (int)nBeams, (long)TEL_SHM_BlockSize, (long)(nBeams * (long)DataSize));
  fprintf(stderr, "-------------------------------------------------------------------------------------\n");

  int done;
  long loop_cntr = 0;
  while (1)
  {
    done = process_acq();
    if (done < 0)
    {
      usleep(100000);
    }
    else if (loop_cntr++ >= obs_duration)
      break;
  }
  return 0;
}
