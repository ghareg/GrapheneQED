CPPFLAGS =

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

jobDetuning: jobDetuning.o basis.o opermask.o chkopts
	-${CLINKER} -o jobDetM1 jobDetuning.o basis.o opermask.o ${SLEPC_EPS_LIB} -lgsl -lgslcblas
	${RM} jobDetuning.o basis.o opermask.o

DATAPATH = ${SLEPC_DIR}/share/slepc/datafiles/matrices

runex:
	-@${MPIEXEC} -n 4 ./jobDetM1 -eps_nev 20 -eps_ncv 50 -eps_tol 1e-4 -eps_smallest_real -eps_type krylovschur > jobDetM1.dat 2>&1;	
