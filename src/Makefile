#------------------------
# basic makefile for mag-ace
#------------------------
#
#  makes:
#  mace
#  libace
#  all
#  clean
#
#--------------------------------------------------------
# adapt settings
#-------------------------------------------------------- 
VER = serial
#VER = parallel

# COMP = intel
COMP = gfortran
# use to add more compilers: else ifeq

BIN = fast

#---------------------------------------------------------
# remainder should not require changes
#----------------------------------------------------------

ifeq ($(BIN),fast)
 PFLAGS2 = NODEBUG
 PFLAGS3 = USELOOKUP
 PFLAGS4 = NOEXTRAVERBOSE
endif

ifeq ($(VER),serial)
 PFLAGS = SERIAL
# INTEL
 ifeq ($(COMP),intel)	
  FC = ifort
  ifeq ($(BIN),fast)
#   CFLAGS = -O2 -xhost
     CFLAGS = -fast
  else
   CFLAGS = 
  endif
 else
# GFORTRAN as default
  FC = gfortran
  ifeq ($(BIN),fast)
   CFLAGS = -Ofast
  else
   # CFLAGS = -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
   CFLAGS = -g -Wuninitialized -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -pedantic -pedantic-errors
  endif
 endif
endif

# parallel
ifeq ($(VER), parallel)
 PFLAGS = MPIPARALLEL
# INTEL
 ifeq ($(COMP),intel)	
  FC = mpif90
  ifeq ($(BIN),fast)
#   CFLAGS = -O2 -xhost
   CFLAGS = -fast
  else
   CFLAGS = 
  endif
 else
# GFORTRAN as default
  FC = mpif90
  ifeq ($(BIN),fast)
   CFLAGS = -Ofast
  else
   # CFLAGS = -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
   CFLAGS = -g -Wuninitialized -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal -pedantic -pedantic-errors
  endif
 endif
endif


LDFLAGS =


# OPEN MP with gfortran
# not to forget: export OMP_NUM_THREADS=8
# CSFLAGS = -fopenmp -O3
# LDFLAGS = -fopenmp

# module load 64/intel/0016.0.047
# FCS = ifort
# CSFLAGS = -fast
# CSFLAGS = -O3
# LDFLAGS = 


#------------------------------------------------------------------
ifeq ($(VER), parallel)
 TARGET = mace.parallel
else
 TARGET = mace.serial
endif
LIB = libace

#-------------------------------------------------------------------

NEILISTDIR = ../neilist
MINSOLVEDIR = ../external
LBFGSDIR = ../external/lbfgs/Lbfgsb.3.0

LAPACKDIR = ../lapack/lapack-3.8.0
LIBLAPACK = ${LAPACKDIR}/liblapack.a ${LAPACKDIR}/librefblas.a
#LIBLAPACK = /usr/lib/x86_64-linux-gnu/lapack/liblapack.a /usr/lib/x86_64-linux-gnu/blas/libblas.a
#LIBLAPACK = -L/usr/lib64 -llapack -lblas

MODS = modules.o

MACE = maceutil.o mace.o

ACE  = ace.o

BASE =  neiutil.o util.o forces.o b_basis.o a_basis.o spherical.o radial.o

NEILIST = ${NEILISTDIR}/neilist.o

MINSOLVE = ${MINSOLVEDIR}/minsolve.o

LBFGS = ${LBFGSDIR}/blas.o ${LBFGSDIR}/linpack.o ${LBFGSDIR}/lbfgsb.o ${LBFGSDIR}/timer.o

OBJS = ${MODS} ${NEILIST} ${MACE} ${BASE} ${MINSOLVE} ${LBFGS}

# LIB needs global NEILIST objs for numerical forces only -> should be disentangled
LIBOBJS = ${MODS} ${NEILIST} ${ACE} ${BASE}

#-----------------------------------------------------------------

DEFAULT: $(TARGET)

all	: $(TARGET) $(LIB)

$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBLAPACK)

$(LIB):	libace.a
libace.a: $(LIBOBJS)
	ar crv $@ $(LIBOBJS) $(LIBLAPACK)
	ranlib libace.a

#------------------------------------------------------------------

# exceptions for external files
${NEILISTDIR}/neilist.o:
	$(FC) -c $(CFLAGS) ${NEILISTDIR}/neilist.f90 -o ${NEILISTDIR}/neilist.o

${MINSOLVEDIR}/minsolve.o:
	$(FC) -c $(CFLAGS) ${MINSOLVEDIR}/minsolve.f -o ${MINSOLVEDIR}/minsolve.o

${LBFGSDIR}/blas.o:
	$(FC) -c $(CFLAGS) ${LBFGSDIR}/blas.f -o ${LBFGSDIR}/blas.o

${LBFGSDIR}/linpack.o:
	$(FC) -c $(CFLAGS) ${LBFGSDIR}/linpack.f -o ${LBFGSDIR}/linpack.o

${LBFGSDIR}/lbfgsb.o:
	$(FC) -c $(CFLAGS) ${LBFGSDIR}/lbfgsb.f -o ${LBFGSDIR}/lbfgsb.o

${LBFGSDIR}/timer.o:
	$(FC) -c $(CFLAGS) ${LBFGSDIR}/timer.f -o ${LBFGSDIR}/timer.o

# regular files and files types
%.o: %.f
	$(FC) -c $(CFLAGS) $<

%.o: %.f90
	$(FC) -c $(CFLAGS) $<

%.o: %.F90
	$(FC) -c $(CFLAGS) -D $(PFLAGS) -D $(PFLAGS2) -D $(PFLAGS3) -D $(PFLAGS4) $<

#--------------------------------------------------------------------------------

clean:
	rm -f *.o *.mod *.MOD ${NEILISTDIR}/*.o ${MINSOLVEDIR}/*.o ${LBFGSDIR}/*.o
	rm -rf *.dSYM
	rm -f *__genmod.f90
	rm -f libace.a

