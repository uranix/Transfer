#CFLAGS=-m32 -O2 -Wall
#CXXFLAGS=$(CFLAGS)

CFLAGS=-m32 -O0 -g -Wall -ImeshProcessor
CXXFLAGS=$(CFLAGS)

CUDA_INSTALL_PATH=/usr/local/cuda
NVCC=$(CUDA_INSTALL_PATH)/bin/nvcc
NVCCFLAGS= -I. -O0 -g -Xcompiler "$(CFLAGS)"\
		   -m32 \
		   -arch sm_20 \
		   --keep --keep-dir cufiles

LDFLAGS=-m32
PTXFLAGS= -v -O2
LIBS= -L/usr/local/cuda/lib/ -LmeshProcessor -lmesh3d -lcudart

OBJS=main.o wrapper.o DeviceAngularData.o DeviceMeshData.o LebedevQuad.o \
	 HemiQuad.o Spherical.o AngularData.o MeshData.o

TARGET=main

.PHONY: all clean deepclean


all: main

%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -c $<

wrapper.o:: kernels.cu

meshProcessor/libmesh3d.a: meshProcessor/*.h meshProcessor/*.h
	$(MAKE) -C meshProcessor all

main: meshProcessor/libmesh3d.a $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

clean: 
	rm -f cufiles/*
	rm -f $(TARGET) $(OBJS)

deepclean: clean
	$(MAKE) -C meshProcessor clean

