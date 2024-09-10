#ifndef SOFTCORR_H
#define SOFTCORR_H

#define DasHeaderKey 2031
#define DasBufferKey 2032

#define DasHeaderKey_SIM 3031
#define DasBufferKey_SIM 3032

#define NCHANNELS 2048
#define NParallel 5
#define NSerial 10
#define NBeams (NParallel * NSerial)
// #define process_psr_SHM_BlockSize (512 * NCHANNELS)
// 32 blocks of the multi-beam version of TEL_SHM (process_psr), each having 800 time samples and 4096 channels; 8-bit integers.
// Time resolution at this stage is 1.31072 ms, which makes one block of the multi-beam version of TEL_SHM (process_psr) 1.048576 s long.
#define FFT_Samps_per_Block 800
#define TEL_SHM_BlockSize (FFT_Samps_per_Block * NCHANNELS)
#define TEL_to_FRB_Block_factor 32                             // 32 blocks of TEL_SHM forms 1 block of FRB_SHM.
#define DataSize (TEL_to_FRB_Block_factor * TEL_SHM_BlockSize) // Size of one block of FRB_BEAM_SHM = 100 MB; 33.554432 s.
#define MaxDataBlocks 12                                       // There will be 12 blocks, each of the above mentioned size, totalling 1.171875 GB; 402.653184 s.
#define Blocks_in_Buf_count 4

typedef struct
{
  unsigned int active, status, is_buf_empty;
  double pc_time, ref_time, rec_time; /***** pc_time = time of day in sec, unix_time = pc_time + ref_time *****/ /***** rec_time = IST *****/
  struct timeval timestamp_gps[MaxDataBlocks];
  double blk_nano[MaxDataBlocks];
  unsigned int flag, curBlock, curRecord, blockSize, nBeams;
  int overFlow;
  unsigned char data[(long)NBeams * (long)(DataSize) * (long)(MaxDataBlocks)];
} DataBuffer;
// DataHeader is the same as TEL_SHM header (BeamHeaderType), the struct definition is given in gmrt_newcorr.h.

#endif // SOFTCORR_H
