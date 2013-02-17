EXEC = gcdp
EXECC = gcdpc
OBJS = gcdp_const.o common.o main.o set_value.o calc_dv_du.o\
       update.o updateu.o stime.o ltime.o output.o starf.o starfd.o\
       cool.o gtree.o tree.o treef.o dmtree.o ldmtree.o dmtreef.o ran.o\
       forcont.o para.o dtyields.o setyield.o\
       atfunc.o esolve.o setcool.o \
       ddecb.o ddecdm.o indexxl.o boundary.o\
       set_valuedm.o\
       setslgtr.o setsldmtr.o kernel.o setmext.o phcurve.o
INIOBJ = set_ini.o
CONTOBJ = cont.o
DEFFILE = gcdp.def
###  for MPICH 
OPTIONS=-O -mcmodel=medium
### kepler 
# MPI_DIR=/cluster/mpich2_3.0
### mac
MPI_DIR = /Users/dkawata/share/mpi
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
FFLAGS = ${OPTIONS} -I${MPI_INC}
LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS} 
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
### for SX-6 isas 20/11/2004
#FORT = sxmpif90
#MPI_INC = 
#MPI_LIB = 
#MLIBS = 
#OPTIONS = -Cvopt -R5 -Wf"-pvctl noassume vwork=stack"
#FFLAGS = ${OPTIONS} 
#LIBS = ${MLIBS}  
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

clean :
	rm *.o  ${EXEC} ${EXECC} ${EXECSC} ${EXEG} ${EXEGC} *.f90 *.mod
