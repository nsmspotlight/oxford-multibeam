CC     = gcc
CFLAGS = -g -Wall -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -fcommon

BEAMSHM_SIM = BeamSHM_Sim.c

OBJ         = simfpsr.o gpu_acqpsr.o  BeamSHM_Sim.o

all : BeamSHM_Sim

BeamSHM_Sim.o :  $(BEAMSHM_SIM)
	$(CC) $(CFLAGS) -c $(BEAMSHM_SIM) -o BeamSHM_Sim.o

gpu_acqpsr.o:gpu_acqpsr.c
	$(CC) $(CFLAGS) -o gpu_acqpsr.o -c gpu_acqpsr.c 

simfpsr.o:simfpsr.c
	$(CC) $(CFLAGS) -o simfpsr.o -c simfpsr.c 

BeamSHM_Sim : $(OBJ)
	$(CC) $(CFLAGS) -o BeamSHM_Sim $(OBJ)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean :
	rm -f *.o 
