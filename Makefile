CUDA_INSTALL_PATH=/usr/local/cuda

CFLAGS=-m32 -O2 -Wall -ImeshProcessor -I"$(CUDA_INSTALL_PATH)/include"
CXXFLAGS=$(CFLAGS)

NVCC=$(CUDA_INSTALL_PATH)/bin/nvcc
NVCCFLAGS= -I. -O2 -Xcompiler "$(CFLAGS)"\
		   -m32 \
		   --keep --keep-dir cufiles

LDFLAGS=-m32 -g
PTXFLAGS= -v -O2
LIBS= -L/usr/local/cuda/lib/ -LmeshProcessor -lmesh3d -lcuda

OBJS=main.o CudaContext.o DeviceAngularData.o DeviceMeshData.o LebedevQuad.o \
	 HemiQuad.o Spherical.o AngularData.o MeshData.o Config.o

CUBIN=kernels.sm_12.cubin kernels.sm_13.cubin kernels.sm_20.cubin

TARGET=main

.PHONY: all clean deepclean


all: main $(CUBIN)

%.sm_12.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -gencode=arch=compute_12,code=sm_12 -cubin $< -o $@

%.sm_13.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -gencode=arch=compute_13,code=sm_13 -cubin $< -o $@

%.sm_20.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -gencode=arch=compute_20,code=sm_20 -cubin $< -o $@

meshProcessor/libmesh3d.a: meshProcessor/*.h meshProcessor/*.h
	$(MAKE) -C meshProcessor all

main: meshProcessor/libmesh3d.a $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

clean: 
	rm -f cufiles/*
	rm -f $(TARGET) $(OBJS) $(CUBIN)

deepclean: clean
	$(MAKE) -C meshProcessor clean

