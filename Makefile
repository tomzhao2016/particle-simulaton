#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
MPCC = CC
OPENMP = -openmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -O3
LIBS =


TARGETS = serial openmp mpi autograder mpi_helper_test

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common2.o 
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common2.o
mpi_helper_test: mpi_helper_test.o common2.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi_helper_test.o common2.o 



autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common2.h mpi_helper.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
common2.o: common2.cpp common2.h
	$(CC) -c $(CFLAGS) common2.cpp
mpi_helper.o: mpi_helper.cpp mpi_helper.h common2.h
	$(CC) -c $(CFLAGS) mpi_helper.cpp
mpi_helper_test.o: mpi_helper_test.cpp mpi_helper.h common2.h
	$(MPCC) -c $(CFLAGS) mpi_helper_test.cpp
clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
