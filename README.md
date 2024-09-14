# `shmutils`

### Shared memory utilities for SPOTLIGHT.

Currently contains the following software:

- `Beam_SHM_SIM`: Emulates the SPOTLIGHT correlator (`TEL_SHM`) by copying data from a file into `TEL_SHM` progressively. This allows one to test software meant for the SPOTLIGHT pipeline on any of their local systems, instead of the actual SPOTLIGHT system itself. Written by [Sanjay Kudale](https://github.com/ksanjay-github).

- `Ares`: Reads beam data from the SPOTLIGHT multi-beam correlator (`TEL_SHM`), and copy it to `FRB_SHM`. Note that this does not process the data in any way. This was built to allow building and testing the rest of the pipeline while a multi-beam version of `gptool` is being built, which is supposed to do the same thing as `Ares`, except that it will mitigate RFI from the data before copying it to `FRB_SHM`. Written by [Kenil Ajudiya](https://github.com/Kenil-Ajudiya).

- `Dump_multi-beam_FRB_SHM`: Reads multi-beam `FRB_SHM`, and dumps the data to disk. Note that you can specify how many seconds of data to dump, and which beam to dump. If the beam number specified is less than 0, all beams are dumped into one file. Written by [Kenil Ajudiya](https://github.com/Kenil-Ajudiya).

Each directory has a `Makefile` that can be used to compile the software, and a `run_*.sh` file that can be used to run it. Note that both might require some editing to make things work.
