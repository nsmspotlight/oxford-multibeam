/*
 *
 */

#ifndef  SOFTCORR_H
#define  SOFTCORR_H

#define  DasHeaderKey   1031
#define  DasBufferKey   1032

#include  <sys/ipc.h>
#include  <sys/shm.h>
#include  <sys/types.h>

#include "gmrt_newcorr.h"

#define DECIMATE 1
#define  BeamBufferSize  (32*1024*1024)/(4*DECIMATE)
#define  VltBeamBufferSize  (256*1024*1024)/(DECIMATE) // to be worked out. 8-bit voltage beam
//#define  VltBeamBufferSize  (128*1024*1024)/(DECIMATE) // to be worked out. 4-bit voltage beam //SSK
#define  IABeamBufferSize  (32*1024*1024)/(4*2*DECIMATE)
#define  MaxDataBlocks  4
//#define  MaxDataBlocks  64
#define  TimeSize       sizeof(double)

#define  DasHeaderKey1   1031
#define  DasBufferKey1   1032

enum BufFlag { _BufMarked=1, _BufReady=1<<1, _BufOverflow=1<< 2, _BufFinish=1<<3, _MaxBufs=100 };

typedef struct
{
  unsigned int flag, curBlock, curRecord, blockSize;
  int overFlow;
  double pcTime[MaxDataBlocks],dataTime[MaxDataBlocks];
  char data[BeamBufferSize*MaxDataBlocks];
} DataBuffer;

typedef struct
{
  unsigned int flag, curBlock, curRecord, blockSize;
  int overFlow;
  double pcTime[MaxDataBlocks],dataTime[MaxDataBlocks];
  short int data[2*IABeamBufferSize*MaxDataBlocks];
} DataBufferIA;

typedef struct
{
  unsigned int active, status;
  char file[100];
  double pcTime, dataTime, refTime;
  struct timeval timestamp[MaxDataBlocks];
  struct timeval timestamp_gps[MaxDataBlocks];
  double blk_nano[MaxDataBlocks];
} DataHeader;

// Voltage beam related.
typedef struct
{
  unsigned int active, status;
  char file[100];
  double pcTime, dataTime, refTime;
  struct timeval timestamp[MaxDataBlocks];
  struct timeval timestamp_gps[MaxDataBlocks];
  double blk_nano[MaxDataBlocks];
  int iteration[MaxDataBlocks];
} DataHeader1; // Voltage-1

typedef struct
{
  unsigned int flag, curBlock, curRecord, blockSize;
  int overFlow;
  double pcTime[MaxDataBlocks],dataTime[MaxDataBlocks];
  char data_p1[(VltBeamBufferSize)*(MaxDataBlocks)]; // 8 MB 
  char data_p2[(VltBeamBufferSize)*(MaxDataBlocks)]; // 8 MB 
  //char data_p1[17179869184]; // 8 MB 
  //char data_p2[17179869184]; // 8 MB 
} DataBuffer1;  // Voltage-1
// Voltage beam related.

#endif  // SOFTCORR_H
