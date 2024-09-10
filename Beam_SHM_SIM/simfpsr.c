/*Renamed previous simfpul.c to simfpsr.c  -- 07/11/1979 */
/* simfpul.c : Dec 6, 1999.
   This program reads data from file and simulates the acquisition of data.

   modified from simf.c of Sep 28, 1998.
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include "acqpsr.h"
#include "pcidev2.h"
#include "parallel_correlator.h"
#include "gsb_unshuf.h"

enum
{
  BlockSize = 16384
};
static int FD = -1, AcqOn = 0, Bufs = 1; //, BlocksRead=0;
static long LastTime[2];
extern double BLKtime;

extern long RefDayTime;
extern double *STAeqn, *STAval;
extern unsigned STA, STAind, *STAseq, STAeqnSet;

#define USec 1000000
#define TDiff(x, y) (((y).tv_sec - (x).tv_sec) * USec + (y).tv_usec - (x).tv_usec)
#define Tdiff(x, y) (((y)[0] - (x)[0]) * USec + (y)[1] - (x)[1])
#define Tadd(t, u)              \
  {                             \
    if (((t)[1] += (u)) > USec) \
    {                           \
      (t)[0] += (t)[1] / USec;  \
      (t)[1] %= USec;           \
    }                           \
  }

int data_open(char *fn, int flag)
{
  if (!fn || !*fn)
    return -1;
  if (FD != -1)
    return -2;
  FD = open(fn, O_RDONLY);
  if (FD < 0)
  {
    perror(fn);
    return -3;
  }
  return FD;
}
/*
int data_close(int fd)
{ extern int set_io(int, unsigned);
  if (fd != FD || FD == -1) return -1;
  if (AcqOn) set_io(FD, 0);
  FD = -1; Bufs = 1;
  return 0;
}
*/
int data_read(int fd, char *buf, size_t size)
{
  int m, n;
  if (fd != FD || FD == -1 || AcqOn == 0)
    return -1;
  if (size < BlockSize)
    return -2;
  n = read(FD, buf, BlockSize);
  if (n < BlockSize)
  {
    lseek(FD, 0, SEEK_SET);
    if (n < 0)
      n = 0;
    m = read(FD, buf + n, BlockSize - n);
    if (m + n < BlockSize)
      return -3;
  }
  BlocksRead++;
  return BlockSize;
}
void compute_time(long *new_time, long add_usec)
{
  /*
  if (random() % 1000 < 10)
    add_usec = add_usec + 2*(random() % add_usec);
  */
  Tadd(new_time, add_usec);
}
int get_block_info(BlockType *blockp)
{
  static int blocks_last = 0;
  //long t_req[2], t_curr[2], usec;
  long t_req[2], usec;
  blockp->status = blockp->block = 0;
  blockp->count = BlocksRead;
  memcpy(t_req, LastTime, 2 * sizeof(long));
  usec = (long)((BlocksRead - blocks_last) * BLKtime * BlockSize / RecSize * 1e6);
  Tadd(LastTime, usec);
  compute_time(t_req, usec);
  blocks_last = BlocksRead;
  /*
    do
    { gettimeofday((void*)t_curr, NULL);
      if (Tdiff(t_req,t_curr) > 0) break;
      usleep(16);
    } while (1);
  */
  memcpy(blockp->time_val, t_req, 2 * sizeof(long));
  return 0;
}
// int set_extra_bufs(int fd, unsigned bufs)
int set_extra_bufs(unsigned bufs)
{
  //  if (fd != FD || FD == -1 || AcqOn == 1) return -1;
  if (bufs > 10)
    return -1;
  Bufs += bufs;
  return 0;
}

// int set_io(int fd, unsigned on)
int set_io(unsigned on)
{
  //  if (fd != FD || FD == -1) return -1;
  if (on)
  {
    if (AcqOn == 1)
      return 0;
    BlocksRead = 0;
    gettimeofday((void *)LastTime, NULL);
    AcqOn = 1;
  }
  else
    AcqOn = 0;
  return 0;
}

// int get_port_info(int fd, PortInfoType *port_info)
int get_port_info(PortInfoType *port_info)
{
  //  if (fd != FD || FD == -1) return -1;
  //  port_info->block_size = BlockSize; /* That is all at the moment */
  port_info->block_size = IABeamBufferSize; /* That is all at the moment */
  return 0;
}

// void print_port_info(int fd, FILE *f)
/*
void print_port_info(FILE *f)
{
//  if (fd != FD || FD == -1) return;
  fprintf(f,"No PortInfo available, BLKtime=%f sec\n", BLKtime);
}
*/
