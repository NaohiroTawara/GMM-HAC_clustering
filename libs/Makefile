# Makefile

CC = g++
CCFLAGS = -O2 
OPENCVLIB = `pkg-config --libs opencv`
OPENCVINC = `pkg-config --cflags opencv`
EIGENINC =  -I$(HOME)/work/spkrCl/libs/eigen-3.2.1
RANDLIB = -L/home/tawara/usr/lib -lrandlibc

LFILE = libspc.a
AR = ar
ARFLAG = crsv


#RANDLIB = -L/home/tawara/usr_64bit/lib -lrandlibc

TARGET = htkdata.o \
         htkmodel.o \
         eval.o \
         evalK.o \
         evalDER.o \
         segment.o

$(LFILE): $(TARGET)
	$(AR) $(ARFLAG) $(LFILE) $(TARGET) 
	@rm -f $(TARGET)

all: $(TARGET)

evalK.o: evalK.cc eval.cc
	$(CC) -c $^ $(CCFLAGS)

evalDER.o: evalDER.cc eval.cc
	$(CC) -c $^ $(CCFLAGS)

eval.o: eval.cc
	$(CC) -c $^ $(CCFLAGS)

htkdata.o: htkdata.cc 
	$(CC) -c $^ $(CCFLAGS)

htkmodel.o: htkmodel.cc 
	$(CC) -c $^ $(CCFLAGS) $(OPENCVINC)

segment.o: segment.cc
	$(CC) -c $^ $(CCFLAGS) $(EIGENINC)

clean:
	rm -rf *.o *.a

