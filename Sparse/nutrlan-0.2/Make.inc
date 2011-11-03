############################################################
# nuTRLan make include file for a plain linux installation with LAPACK,
# BLAS and MPI in /usr/local
############################################################
#
# C compiler
  CC      = cc
  MPICC   = mpicc
  LOADER  = cc
  PLOADER = mpicc
#
# Debuging
  DEBUG = -Wall
  OPT = -O3 -DBLAS -ggdb
#
# Make command
  MAKE = make
#  MAKEFLAGS = TOP=${PWD}
#
# nuTRLan home directory
#  TOP    = ${PWD}
  NUTRLAN = $(TOP)/libnutrlan.a
#
# External directories
  DIR_CLAPACK = 
#
# External libraries
  LIB_DIR =
  MATH    =
  BLAS    = -lblas
  LAPACK  = -llapack
  CLAPACK = 
  F2C     = 
  LIB     = $(LIB_DIR) $(LAPACK) $(BLAS) $(CLAPACK) $(F2C) $(MATH) -lpthread
  PLIB    = $(LIB) -lmpi 
#
# Include files
  INC_CLAPACK = 
  INC_NUTRLAN  = -I$(TOP)/INC
  INC         = $(INC_CLAPACK) $(INC_NUTRLAN)
#
# Object files
  OBJ_REAL = dsort2.o dstqrb.o trlan.o trlaux.o trlcore.o restart.o trlmap.o
  OBJ_CPLX = zdgemm.o zdgemv.o zdaxpy.o ztrlan.o ztrlaux.o ztrlcore.o
  OBJ = $(OBJ_REAL) $(OBJ_CPLX)
  OBJ_SERIAL = $(OBJ) trl_comm_none.o ztrl_comm_none.o
  OBJ_PARALLEL = $(OBJ) trl_comm_mpi.o ztrl_comm_mpi.o
