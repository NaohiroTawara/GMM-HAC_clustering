# Makefile                                                                      

CC = g++
CCFLAGS = -O2
MPFLAGS = -fopenmp
#DEBUG = -pg
MPFLAGS = -fopenmp
OPENCVINC = `pkg-config --cflags opencv`
OPENCVLIB = `pkg-config --libs opencv`
#EIGENINC =  -I/usr/local/include/eigen3/
EIGENINC =  -I$(HOME)/work/spkrCl/libs/eigen-3.2.1/

SCDIR = $(HOME)/work/spkrCl/libs/tags/1.0.5
#SCDIR = $(HOME)/work/spkrCl/libs/tags/1.0.5.1 # for old version of diag
SCINC = -I$(SCDIR)
SCLIB = -L$(SCDIR) -lspc 
SCOBJ = gmm.o gmmcluster.o gmmclustering.o spkr_cl_GMM-HAC

TARGET = spkr_cl_GMM-HAC gmmmlcluster.o gmmclustering.o gmm.o

all: $(TARGET)

spkrLib:
	make -C $(SCDIR)

spkr_cl_GMM-HAC: gmm.o gmmmlcluster.o gmmclustering.o main.cc 
	$(CC) -o $@ $^  $(CCFLAGS) $(EIGENINC) $(SCINC) $(SCLIB) $(MPFLAGS)

gmmmlcluster.o: gmmmlcluster.cc
	$(CC) -c $^ $(CCFLAGS) $(EIGENINC) $(SCINC)

gmmclustering.o: gmmclustering.cc
	$(CC) -c $^ $(CCFLAGS) $(EIGENINC) $(SCINC)

gmm.o: gmm.cc gmm.h
	$(CC) -c $^ $(CCFLAGS) $(EIGENINC) $(SCINC) $(MPFLAGS)

clean:
	rm -rf *.o

#end of file  

