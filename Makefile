EXEC = gcdp
EXECC = gcdpc
OBJS = gcdp_const.o common.o main.o set_value.o calc_dv_du.o\
       update.o updateu.o stime.o ltime.o output.o starf.o starfd.o\
       cool.o gtree.o tree.o treef.o dmtree.o ldmtree.o dmtreef.o ran.o\
       forcont.o para.o dtyields.o setyield.o\
       atfunc.o esolve.o setcool.o \
       ddecb.o ddecdm.o indexxl.o boundary.o\
       set_valuedm.o set_chi.o\
       setsldmtr.o kernel.o setmext.o allocate.o\
       mesh_gen.o mesh_setrho.o fftw3_init.o mesh_fftf.o mesh_setfp.o erfcc.o\
       phcurve.o phcurven.o mesh_isogpot.o del_dmpart.o set_testpdm.o
INIOBJ = set_ini.o
CONTOBJ = cont.o
DEFFILE = gcdp.def gcdp_const.F95 common.F95
###  for MPICH 
OPTIONS=-g -mcmodel=medium
### kepler 
#MPI_DIR=/cluster/mpich2_3.0
#FFTW_INC=/khome/dka/share/include
#FFTW_LIB=/khome/dka/share/lib
### mac
MPI_DIR = /usr/local
FFTW_INC=/usr/local/include
FFTW_LIB=/usr/local/lib
### legion@UCL
# module unload compilers/intel/2015/update2 
# module unload mpi/intel/2015/update3/intel
# module load compilers/gnu/4.9.2
# module load mpi/openmpi/1.8.4/gnu-4.9.2
# module load fftw/3.3.4-ompi/gnu-4.9.2 
#
# MPI_DIR=/shared/ucl/apps/openmpi/1.8.4/gnu-4.9.2
# OPTIONS=-O3
# FFTW_INC=/shared/ucl/apps/fftw/3.3.4-ompi/gnu-4.9.2/include
# FFTW_LIB=/shared/ucl/apps/fftw/3.3.4-ompi/gnu-4.9.2/lib
CPP = cpp
FORT = ${MPI_DIR}/bin/mpif90
MPI_INC = ${MPI_DIR}/include
MPI_LIB = ${MPI_DIR}/lib 
FFLAGS = ${OPTIONS} -I${MPI_INC} -I${FFTW_INC}
LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS} -L${FFTW_LIB} -lfftw3_mpi -lfftw3 
### for XT4
# FORT=ftn
# CPP=cpp
# OPTIONS= -fastsse -O4 -Mvect=sse
# FFLAGS = ${OPTIONS}
# LIBS=${LINK_MPI}
### cosmos
# FORT=ifort
# CPP=cpp
# OPTIONS=-g -O3 -xHost -align -ansi-alias -mcmodel=medium -traceback
# FFLAGS = ${OPTIONS}
# LIBS=-lmpi -lfftw3_mpi -lfftw3
### darwin
# module add fftw/intel/3.3.3
# FORT=mpif90
# CPP=cpp
# FFLAGS= -O3 -ipo -xHost -lfftw3 -lfftw3_mpi -lfftw3
#
### leicester
#
# module load intel/compilers/14.0.1
# module load intel/impi/4.1.1
#
#OPTIONS=-xAVX -O3 -mcmodel=medium -assume buffered_io -ipo
#OPTIONS=-xAVX -O3 -mcmodel=medium -assume buffered_io
#FFTW_INC=/cm/shared/apps/fftw/intel/3.3.4/include
#FFTW_LIB=/cm/shared/apps/fftw/intel/3.3.4/lib
#MPI_DIR=/cm/shared/apps/ics/impi/4.1.1.036//intel64/bin
#FORT=mpiifort
#CPP=cpp
#LDFLAGS=-Wl,-R,/cm/shared/apps/ics/composer_xe_2013_sp1.1.106/compiler/lib/intel64
#OPTIONS=-O3 -mcmodel=medium
#FFLAGS = ${OPTIONS} -I${FFTW_INC}
#LIBS= ${OPTIONS} -L${FFTW_LIB} -lfftw3_mpi -lfftw3 ${LDFLAGS} 
#
### iridis
# module load openmpi/1.3.3/intel
# module load fftw/3.3.1/intel/double
#FORT=mpif90
#CPP=cpp
#FFTW_INC=/local/software/rh53/fftw/3.3.1/intel/double/include
#FFTW_LIB=/local/software/rh53/fftw/3.3.1/intel/double/lib
#OPTIONS= -O3
#FFLAGS = ${OPTIONS} -I${FFTW_INC}
#LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS} -L${FFTW_LIB} -lfftw3_mpi -lfftw3
#
### guillimin
# module add ifort_icc/14.0.1 mvapich2/1.9-intel FFTW/3.3-mvapich2-intel
# 
# FORT=mpif90
# CPP=cpp
# OPTIONS=-O3 -xHost -align -ansi-alias -mcmodel=medium -traceback
# FFLAGS = -I${MPI_INC} -I${FFTW_INC} ${OPTIONS}
# MPI_DIR=/software/CentOS-6/tools/mvapich2-1.9-intel
# MPI_INC = ${MPI_DIR}/include
# MPI_LIB = ${MPI_DIR}/lib
# FFTW_INC=/software/CentOS-6/libraries/FFTW-3.3/mvapich2-intel/include
# FFTW_LIB=/software/CentOS-6/libraries/FFTW-3.3/mvapich2-intel/lib
# LIBS = -L${MPI_LIB} ${OPTIONS} -L${FFTW_LIB} -lfftw3_mpi -lfftw3
#
# CINECA/EURORA
#
# module load intel/cs-xe-2013--binary openmpi/1.6.5--intel--cs-xe-2013--binary fftw/3.3.3--openmpi--1.6.5--intel--cs-xe-2013--binary
#
#FORT=mpif90
#CPP=cpp
#OPTIONS= -O3
#FFLAGS = ${OPTIONS} -I${FFTW_INC}
#LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS} -L${FFTW_LIB} -lfftw3_mpi -lfftw3
#
# CINECA/PLX
# module load intel/cs-xe-2013--binary intelmpi
#
# FORT=mpif90
# CPP=cpp
# FFTW_INC=/cineca/prod/libraries/fftw/3.3.3/intelmpi--4.1.2--binary/include
# FFTW_LIB=/cineca/prod/libraries/fftw/3.3.3/intelmpi--4.1.2--binary/lib
# OPTIONS= -O3
# FFLAGS = ${OPTIONS} -I${FFTW_INC}
# LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS} -L${FFTW_LIB} -lfftw3_mpi -lfftw3

.SUFFIXES : .o .F95 .f90

${EXEC} : ${OBJS} ${INIOBJ}
	${FORT} -o ${EXEC} ${OBJS} ${INIOBJ} ${FFLAGS} ${LIBS}
	mv ${EXEC} ..

${EXECC} : ${OBJS} ${CONTOBJ}
	${FORT} -o ${EXECC} ${OBJS} ${CONTOBJ} ${FFLAGS} ${LIBS}
	mv ${EXECC} ..

.F95.o: gcdp_const.f90 common.f90
# for intel fortran and XT4
#	${FORT} ${CPPFLAGS} ${FFLAGS} -c $<
# for gnu
	${CPP} -P $< > $*.f90
	${FORT} ${FFLAGS} -c $*.f90

${OBJS} : ${DEFFILE}
${INIOBJ}: ${DEFFILE}
${CONTOBJ} : ${DEFFILE}
${CONTSOBJ} : ${DEFFILE}
%.o : %.mod

clean :
	rm *.o  ${EXEC} ${EXECC} ${EXECSC} ${EXEG} ${EXEGC} *.f90 *.mod
