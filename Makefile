CUDA_INSTALL_PATH=/usr/local/cuda
NVCC=$(CUDA_INSTALL_PATH)/bin/nvcc
NVCCFLAGS= -I. -O2 \
		   -m32 \
		   -arch sm_13 \
		   --keep --keep-dir cufiles

#CFLAGS=-m32 -O2 -Wall
#CXXFLAGS=$(CFLAGS)

CFLAGS=-m32 -O0 -g -Wall -ImeshProcessor
CXXFLAGS=$(CFLAGS)

LDFLAGS=-m32
PTXFLAGS= -v -O2
LIBS= -L/usr/local/cuda/lib/ -LmeshProcessor -lcudart -lmesh3d

OBJS=main.o wrapper.o DeviceAngularData.o DeviceMeshData.o LebedevQuad.o \
	 HemiQuad.o Spherical.o AngularData.o MeshData.o
TARGET=main

.PHONY: all clean deepclean

all: main

%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -c $<

meshProcessor/libmesh3d.a:
	$(MAKE) -C meshProcessor all

main: meshProcessor/libmesh3d.a $(OBJS)
	gcc $(LDFLAGS) -o $@ $^ $(LIBS)

clean: 
	rm -f cufiles/*
	rm -f $(TARGET) $(OBJS)

deepclean: clean
	$(MAKE) -C meshProcessor clean

