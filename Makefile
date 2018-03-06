# Load CUDA using the following command
# module load cuda
#
CC = nvcc
CFLAGS = -O3 -arch=compute_37 -code=sm_37
NVCCFLAGS = -O3 -arch=compute_37 -code=sm_37
LIBS = 

TARGETS = serial gpu autograder

all:	$(TARGETS)

serial: serial.o common2.o
	$(CC) -o $@ $(LIBS) serial.o common2.o
gpu: gpu.o common2.o
	$(CC) -o $@ $(NVCCLIBS) gpu.o common2.o
autograder: autograder.o common2.o
	$(CC) -o $@ $(LIBS) autograder.o common2.o

serial.o: serial.cu common2.h
	$(CC) -c $(CFLAGS) serial.cu
autograder.o: autograder.cu common2.h
	$(CC) -c $(CFLAGS) autograder.cu
gpu.o: gpu.cu common2.h
	$(CC) -c $(NVCCFLAGS) gpu.cu
common2.o: common2.cu common2.h
	$(CC) -c $(CFLAGS) common2.cu

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
