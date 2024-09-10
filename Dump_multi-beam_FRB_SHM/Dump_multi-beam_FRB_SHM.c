/******************************************************************************************************************************************************
 * This program reads the multi-beam raw data from the multi-beam FRB_SHM of SPOTLIGHT and dumps it into a user-specified .raw file.
 * COMPILE: gcc -O3 -Wall -Wextra -std=gnu99 Dump_multi-beam_FRB_SHM.c -o Dump_multi-beam_FRB_SHM
 * USAGE: ./Dump_multi-beam_FRB_SHM -f [file path to dump raw data] -d [duration (in seconds)] {-b [index of the beam to dump]} {-h}
 * Author: Kenil Ajudiya. Contact: kenilr@iisc.ac.in                                    Date: Sept 08, 2024
 ******************************************************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/shm.h>
#include <sys/time.h>
#include "frb_shm.h"

DataBuffer *dataBuffer_FRB = NULL;
unsigned int RecNum = 0, DataSeq = 0, nBeams = 0;

int initialize_FRB_SHM()
{
    int id_Data_buffer = shmget(DasBufferKey, sizeof(DataBuffer), SHM_RDONLY);
    if (id_Data_buffer == -1)
    {
        perror("shmget() error");
        return -1;
    }

    dataBuffer_FRB = (DataBuffer *)shmat(id_Data_buffer, 0, SHM_RDONLY);
    if ((void *)dataBuffer_FRB == (void *)-1)
    {
        perror("shmat() error");
        return -2;
    }

    RecNum = dataBuffer_FRB->curRecord;
    DataSeq = dataBuffer_FRB->curBlock;
    nBeams = dataBuffer_FRB->nBeams;
    fprintf(stderr, "nBeams in the FRB_SHM: %u\n", nBeams);
    return 0;
}

void help()
{
    fprintf(stderr, "# HELP: This program reads the multi-beam raw data from the multi-beam FRB_SHM of SPOTLIGHT and dumps it into a user-specified .raw file.\n");
    fprintf(stderr, "# USAGE: ./Dump_multi-beam_FRB_SHM -f [file path to dump raw data] -d [duration (in seconds)] {-b [index of the beam to dump]} {-h}\n");
    return;
}

int main(int argc, char *argv[])
{
    int opt = 0, blk = 1, nBlocks = 0, beam_to_dump = -1;
    FILE *file = NULL;

    while ((opt = getopt(argc, argv, "f:d:b:h")) != -1)
    {
        switch (opt)
        {
        case 'f':
            file = fopen(optarg, "wb");
            if (file == NULL)
            {
                perror("# ERROR: Failed to open the output data file");
                return -1;
            }
            break;
        case 'd':
            nBlocks = atoi(optarg) / (TEL_to_FRB_Block_factor * FFT_Samps_per_Block * 1.31072 / 1000); // nBlocks = obs. dur. / time per block (33.554432 sec)
            if (nBlocks <= 0)
            {
                fprintf(stderr, "# ERROR: Observation duration (in seconds) must be positive. Exiting.\n");
                return -1;
            }
            break;
        case 'b':
            beam_to_dump = atoi(optarg);
            if (beam_to_dump < 0)
                fprintf(stderr, "# WARNING: Index of the beam to be dumped must be positive. Dumping all the beams.\n");
            break;
        default:
            help();
            exit(1);
            break;
        }
    }

    if (initialize_FRB_SHM() < 0)
    {
        fprintf(stderr, "# ERROR: Could not initialize the FRB_SHM to read buffer.\n");
        fprintf(stderr, "# DEBUG HELP: Check if the frb_shm.h file is up-to-date.\n");
        return -1;
    }
    else
        fprintf(stderr, "# LOG: FRB_SHM initialization successful.\n");

    if (beam_to_dump < 0) // dump all the beams
    {
        while (blk <= nBlocks)
        {
            fprintf(stderr, "# LOG: Reading block %d / %d.\n", blk, nBlocks);
            while (DataSeq >= dataBuffer_FRB->curBlock)
            {
                fprintf(stderr, "# LOG: Sleeping peacefully...\n");
                sleep(2);
            }
            // Dumping a record of (NBeams * DataSize) bytes from the buffer to the output file.
            fwrite(dataBuffer_FRB->data + RecNum * (long)DataSize * (long)NBeams, sizeof(unsigned char), DataSize * nBeams, file);

            RecNum = (RecNum + 1) % MaxDataBlocks;
            DataSeq++;
            blk++;
        }
    }
    else // dump only the specified beam
    {
        while (blk <= nBlocks)
        {
            fprintf(stderr, "# LOG: Reading block %d / %d.\n", blk, nBlocks);
            while (DataSeq >= dataBuffer_FRB->curBlock)
            {
                fprintf(stderr, "# LOG: Sleeping peacefully...\n");
                sleep(2);
            }
            // Dumping a record of (NBeams * DataSize) bytes from the buffer to the output file.
            fwrite(dataBuffer_FRB->data + (long)RecNum * (long)DataSize * (long)NBeams + (long)DataSize * (long)beam_to_dump, sizeof(unsigned char), DataSize, file);

            RecNum = (RecNum + 1) % MaxDataBlocks;
            DataSeq++;
            blk++;
        }
    }
    fprintf(stderr, "# LOG: Done reading all the blocks. Closing the output file...\n");

    fclose(file);
    return 0;
}
