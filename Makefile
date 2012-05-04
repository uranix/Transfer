CUDA_INSTALL_PATH=/usr/local/cuda

CFLAGS=-m32 -O0 -g -Wall -ImeshProcessor -I"$(CUDA_INSTALL_PATH)/include"
CXXFLAGS=$(CFLAGS)

NVCC=$(CUDA_INSTALL_PATH)/bin/nvcc
NVCCFLAGS= -I. -O2 -Xcompiler "$(CFLAGS)"\
		   -m32 \
		   -arch sm_13 \
		   --keep --keep-dir cufiles

LDFLAGS=-m32 -g
PTXFLAGS= -v -O2
LIBS= -L/usr/local/cuda/lib/ -LmeshProcessor -lmesh3d -lcuda

OBJS=main.o CudaContext.o DeviceAngularData.o DeviceMeshData.o LebedevQuad.o \
	 HemiQuad.o Spherical.o AngularData.o MeshData.o Config.o

CUBIN=kernel.cubin

TARGET=main

.PHONY: all clean deepclean


all: main $(CUBIN)

%.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -cubin $<

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

