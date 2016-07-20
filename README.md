# CALPUFF_with_OMP
**Developing CALPUFF(v7.2) with OMP so that it can run parallel.**

# Brief description
Base on **_CALPUFF v7.2_**, codes were transformed from fixed format(.for) to free format(.f90).

Parallel calculating tasks are held among PUFF(s), thus PARALLEL region was created at *puff loop* in `subroutine comp` in the file `calpuff.f90`.

Declarations of variables involved in the functions or subroutines invoked in the PARALLEL region, and those in `subroutine comp` were added, so that some OMP attributes, eg. THREADPRIVATE, FISRTPRIVATE, etc., can be defined.

All modifications by me were marked with *wangzhm* in the code files.

Then the codes were compiled with *Intel Visual Fortran(IVF)* with *OMP* on.
The recompiling command and the options related can be found in `CPL_CALPUFF.BAT` and `cpl_unix.bat`.

For more details, please contact me via E-mail: *dreaming@ixy.info*

**Welcome aboard CALPUFFers!**


# License
This repository uses the *Apache License 2.0*.
