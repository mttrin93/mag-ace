# mag-ace

`mag-ace` is the fortran implementation of magnetic atomic cluster expansion.
It provides the basis functionality for running a magnetic ACE potential .

If you are using this code, please cite:

M. Rinaldi, M. Mrovec, A. Bochkarev, Y. Lysogorskiy and R. Drautz, Non-collinear Magnetic Atomic Cluster Expansion for Iron, arXiv preprint arXiv:2305.15137

R. Drautz, Atomic cluster expansion of scalar, vectorial, and tensorial properties including magnetism and charge transfer, Physical Review B 102 (2), 024104 (2020)

### Installation

1. cp `pair_ace.cpp` and `pair_ace.h` in your `LAMMPS/src/SPIN` directory
2. in your ACE directory do `make libace` to generate `libace.a`
3. in `LAMMPS/src/Makefile` add `PathToLibace/libace.a` after `OBJ = ...`
4. in `LAMMPS/src/MAKE` edit `Makefile.serial` or `Makefile.mpi`. After `LINKFLAGS = ...` add `-lstdc++ -lgfortran`, possibly drop `-cxxlib`
5. compile LAMMPS with the `SPIN` package
