EXEC = gcdp
EXECC = gcdpc
OBJS = gcdp_const.o common.o main.o set_value.o calc_dv_du.o\
       update.o updateu.o stime.o ltime.o output.o starf.o starfd.o\
       cool.o gtree.o tree.o treef.o dmtree.o ldmtree.o dmtreef.o ran.o\
       forcont.o para.o dtyields.o setyield.o\
       atfunc.o esolve.o setcool.o \
       ddecb.o ddecdm.o indexxl.o boundary.o\
       set_valuedm.o set_chi.o\
       setslgtr.o setsldmtr.o kernel.o setmext.o phcurve.o allocate.o\
       mesh_gen.o mesh_setrho.o fftw3_init.o mesh_fftf.o mesh_setfp.o erfcc.o
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
MPI_DIR = /Users/dkawata/share/mpi
FFTW_INC=/usr/local/include
FFTW_LIB=/usr/local/lib
### legion@UCL
# MPI_DIR=/usr/mpi/qlogic
# OPTIONS=-O3 -xSSSE3
### guillimin
# MPI_DIR=/software/tools/openmpi-1.4.3-intel/
# OPTIONS=-O3 -ldl
#
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
# LIBS=-lmpi
### darwin
# module add fftw/intel/3.3.3
# FORT=mpif90
# CPP=cpp
# FFLAGS= -O3 -ipo -xHost -lfftw3
#
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
