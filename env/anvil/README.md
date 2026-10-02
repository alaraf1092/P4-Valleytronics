# P4 Anvil + VASP Environment

## System
Cluster: Purdue Anvil
Slurm partition: shared
Slurm account: nnt250008
Anvil username: x-aaraf

## VASP
VASP version: 6.3.2

## Toolchain
GCC: 11.2.0
OpenMPI: 4.1.6
Intel MKL: 2019.5.281

## CPU target
VASP_TARGET_CPU: -march=znver3

## Build
make DEPS=1 -j8 all

## Executables
vasp_std
vasp_gam
vasp_ncl

## Runtime
Before launching VASP, load:
- gcc/11.2.0
- openmpi/4.1.6
- intel-mkl/2019.5.281

## Notes
The working VASP build was compiled through Slurm on Anvil.
The successful Slurm build completed with ExitCode 0.
