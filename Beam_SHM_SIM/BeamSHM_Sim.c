/**********************************************************************
 * Software to emulate real-time multi-beam SHM exactly similar to
 * SPOTLIGHT correlator's multi-beam SHM.
 * USAGE :  ./bin/BeamSHM_Sim <data_file> <nBeams> <nChans>
 *
 * It reads one buffer of 1.048576 sec at a time and put in SHM.
 * Next buffer read exactly after 1.048576 sec, hence real-time effect.
 * Data file once finished, rewinds back to begining.
 *                                        Sanjay Kudale
 *                                        29 May 2024
 * compile : Makefile
 * compile : gcc -o BeamSHM_Sim BeamSHM_Sim.c -D_LARGEFILE64_SOURCE  -D_FILE_OFFSET_BITS=64
 * ********************************************************************/
// #define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include "parallel_correlator.h"
#include <string.h>

#define  BeamHdrKey      1050

int fill_header (BeamHeaderType *beamHeader) {
	FILE *ra_dec_BMSTEER;
  float bm_RA, bm_DEC;
  int  bm_index, bm_subindex;
	int steered_BM_ID;

	beamHeader->BeamGenHdr.BeamHostID= 0;
	strcpy((beamHeader->BeamGenHdr.BeamHostName),"gpbcorr1");
	beamHeader->BeamGenHdr.BeamType[0] = 2; // Beam-type =2 is post-correlation (PC)
	beamHeader->BeamGenHdr.NStokes[0]  = 1; // NStokes = 1 is single stokes, total intensity
	beamHeader->BeamGenHdr.OutputDataFormat=8; // 8-bit data.
	strcpy((beamHeader->ScanTab[0].source.object),"B1929+10");
	beamHeader->ScanTab[0].source.ra_app = 5.119987;  // RA in radians
	beamHeader->ScanTab[0].source.dec_app = 0.1928;  // RA in radians
	beamHeader->ScanTab[0].source.freq[0] = 550.0;  // Frequency of first channel
	beamHeader->corr.corrpar.clock = 400000000 ; // 400 MHz clock
	beamHeader->corr.corrpar.channels = 4096; // 4k channels
  beamHeader->ScanTab[0].source.net_sign[0] = -1; // Band 4, USB data
  // beamHeader->ScanTab[0].source.net_sign[0] = 1; // Band 3, LSB data
	beamHeader->corr.corrpar.f_step = (-1*(beamHeader->ScanTab[0].source.net_sign[0]))*beamHeader->corr.corrpar.clock/(beamHeader->corr.corrpar.channels * 2) ;  //Hz.
	beamHeader->corr.daspar.gsb_acq_bw = 200.0; // 200 MHz Bandwidth
	beamHeader->corr.daspar.gsb_final_bw = 1; // 200 MHz Bandwidth
  beamHeader->BeamGenHdr.PostTimeInt[0]=1;
  beamHeader->BeamGenHdr.PostFreqInt[0]=1;
	beamHeader->BeamGenHdr.SampInterval=1310.72; // micro-second
	//beamHeader->BeamGenHdr.SampInterval= corr.daspar.gsb_final_bw*BeamGen.SampInterval*BeamGen.PostTimeInt[BeamID]/(corr.corrpar.clock/1e6);

	// 2000 beam RA/DEC etc.
	beamHeader->BeamGenHdr.BeamSteeringParams.nPCBaselines=6; // Number of baselines inlcuded in PC
	beamHeader->BeamGenHdr.BeamSteeringParams.nBeamHosts=2; // Number of baselines inlcuded in PC
	beamHeader->BeamGenHdr.BeamSteeringParams.nSteeringBeams=100; // Number of baselines inlcuded in PC
	beamHeader->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode=50; // Number of baselines inlcuded in PC

	if((ra_dec_BMSTEER=fopen("./ra-dec_centers.txt", "r")) == NULL) {
		perror("BM STEER RA DEC file opening problem");
		return (-1);
	} else {
		for ( steered_BM_ID=0; steered_BM_ID<beamHeader->BeamGenHdr.BeamSteeringParams.nSteeringBeams;steered_BM_ID++) {
			fscanf(ra_dec_BMSTEER, "%f %f %d %d\n", &bm_DEC, &bm_RA, &bm_index, &bm_subindex);

			beamHeader->BeamGenHdr.BeamSteeringParams.RA[steered_BM_ID]=bm_RA;
			beamHeader->BeamGenHdr.BeamSteeringParams.DEC[steered_BM_ID]=bm_DEC;
			beamHeader->BeamGenHdr.BeamSteeringParams.Beam_index[steered_BM_ID]=bm_index;
			beamHeader->BeamGenHdr.BeamSteeringParams.Beam_subindex[steered_BM_ID]=bm_subindex;
		}
		fclose(ra_dec_BMSTEER);
	}

	return 0;
}

int main(int argc, char *argv[])
{
  int beam_ShmId;
  int nstokes_beam = 1;
  int beamid = 0;
  int beam_int = 64;
  int samples_per_word = 1;
  int beam_bits = 2;
  int num_beams;
  int two_stage_size = 4;
  int run_count = 0;
  int infile_fd;
  char *indata;
  long nbeam_blocksize = 163840000;
  struct timeval timestamp;
  int nano_sec = 0;
  int add_usec = 48576;
  long nread;
	int nfft;
	int nchan;
	int nbeamsPerHost;

  unsigned long long lasttime, now;
  unsigned long long nano = 1000000000, wait = 1048576;
  struct timespec tm;
  int sleep_tm;
	static int shmBeamHdr;
  BeamHeaderType *beamHeader;

  if (argc != 4)
  {
    fprintf(stderr, "USAGE :  %s <data_file> <nBeams> <nChans>\n", argv[0]);
    return -1;
  }

	nbeamsPerHost = atoi(argv[2]);
	nchan = atoi(argv[3]);
	nfft = 800;
  num_beams = 2 * nbeamsPerHost;
	nbeam_blocksize = nchan * nfft * nbeamsPerHost;

  indata = malloc(nbeam_blocksize);

  shmBeamHdr = shmget (BeamHdrKey, sizeof (BeamHeaderType), 0644 | IPC_CREAT);
  if (shmBeamHdr < 0) { perror ("SHMGET BEAM HEADER"); return -1; }

  beamHeader = (BeamHeaderType *) shmat (shmBeamHdr, 0, 0);
  if (beamHeader == (void *) -1) { perror ("BeamHeader:SHMAT"); return -1; }

  fill_header (beamHeader);
  
  beamHeader->BeamGenHdr.BeamSteeringParams.nSteeringBeamsPerNode = nbeamsPerHost;

  if ((beam_ShmId = initialize_IAPA(nstokes_beam, nchan, beam_int / samples_per_word, beam_bits, beamid, num_beams / (two_stage_size / 2))) < 0)
  {
    perror("IAPA SHMGET ERROR:");
    exit(-1);
  }

  if ((infile_fd = open(argv[1], O_RDONLY | O_LARGEFILE)) == -1)
  {
    perror(argv[1]);
    return -1;
  }
  gettimeofday(&timestamp, 0);
  while (1)
  {
    clock_gettime(CLOCK_REALTIME, &tm);
    lasttime = (tm.tv_nsec + tm.tv_sec * nano);

    if ((nread = read(infile_fd, indata, nbeam_blocksize)) != nbeam_blocksize)
    {
      perror("REWINDING TO BEGINING");
      lseek(infile_fd, SEEK_SET, 0);
      if ((nread = read(infile_fd, indata, nbeam_blocksize)) != nbeam_blocksize)
      {
        perror("READ ERROR on Rewind");
      }
    }

    if ((WriteBeam2SHM(1, timestamp, timestamp, nano_sec, indata, run_count, nbeam_blocksize, num_beams / (two_stage_size / 2)) < 0))
    {
      perror("BEAM WRITING ERROR");
      exit(-1);
    }
    /*

       local_t = localtime((const time_t *) &timestamp);
       frac_sec = timestamp.tv_usec;
       fprintf(stderr, "# %02d:%02d:%02d.%06d%03d Iteration = %d\n", local_t->tm_hour, local_t->tm_min, local_t->tm_sec, (int)frac_sec, (int)nano_sec, run_count);
    */

    clock_gettime(CLOCK_REALTIME, &tm);
    now = tm.tv_nsec + tm.tv_sec * nano;

    timestamp.tv_usec += add_usec;
    if (timestamp.tv_usec >= 1000000)
    {
      timestamp.tv_sec += 2;
      timestamp.tv_usec -= 1000000;
    }
    else
      timestamp.tv_sec += 1;

    run_count++;

    sleep_tm = lasttime / 1000 + wait - now / 1000;
    usleep(sleep_tm);
  }
}
