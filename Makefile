# Touch util.h if you want only to change REAL type
ifeq ($(USE_DOUBLE),1)
    CFLAGS=-DUSE_DOUBLE 
endif

CFLAGS+= -m32 -O2 -Wall -ImeshProcessor

ifeq ($(OS),Windows_NT)
    CUDA_INSTALL_PATH=$(CUDA_PATH)
    CCBIN=-ccbin "D:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin"
    CFLAGS+= -I"$(CUDA_INSTALL_PATH)include\\"
    CCBINFLAGS+= -W3 -ImeshProcessor
    LIBS= -L"$(CUDA_INSTALL_PATH)lib\Win32" -LmeshProcessor -lmesh3d -lcuda
    NVCC="$(CUDA_INSTALL_PATH)bin\nvcc.exe"
    TARGET=main.exe
else
    CUDA_INSTALL_PATH=/usr/local/cuda
    CCBIN=
    CFLAGS+= -I"$(CUDA_INSTALL_PATH)/include"
    CCBINFLAGS+=$(CFLAGS)
    LIBS= -L"$(CUDA_INSTALL_PATH)/lib/" -LmeshProcessor -lmesh3d -lcuda
    NVCC=$(CUDA_INSTALL_PATH)/bin/nvcc
    TARGET=main
endif

CXXFLAGS=$(CFLAGS)

NVCCFLAGS= $(CCBIN) -I. -O2 -Xcompiler "$(CCBINFLAGS)"\
		   -m32 \
		   --keep --keep-dir cufiles
LDFLAGS=-m32 -g
PTXFLAGS= -v -O2

OBJS=main.o CudaContext.o DeviceAngularData.o DeviceMeshData.o LebedevQuad.o \
	 HemiQuad.o Spherical.o AngularData.o MeshData.o Config.o Util.o Eigenvalue.o \
	 E1.o

DEPS=$(OBJS:.o=.dep) kernels.dep

CUBIN=kernels.sm_12.cubin kernels.sm_13.cubin kernels.sm_20.cubin

.PHONY: all clean deepclean cubin libs cleandeps 

all: libs cubin $(TARGET)

Makefile.deps: cleandeps $(DEPS)

cleandeps: 
	rm -f Makefile.deps

%.dep : %.cpp
	g++ $(CXXFLAGS) -M $< >> Makefile.deps

%.dep : %.cu
	$(NVCC) $(NVCCFLAGS) -M $< | sed "s/kernels.o/$(CUBIN)/">> Makefile.deps

cubin: $(CUBIN)

libs: 
	$(MAKE) -C meshProcessor all

%.sm_12.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -gencode=arch=compute_12,code=sm_12 -cubin $< -o $@

%.sm_13.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -gencode=arch=compute_13,code=sm_13 -cubin $< -o $@

%.sm_20.cubin : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -gencode=arch=compute_20,code=sm_20 -cubin $< -o $@

$(TARGET): $(OBJS)
	g++ $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean: 
	rm -f cufiles/*
	rm -f $(TARGET) $(OBJS) $(CUBIN)

deepclean: clean
	$(MAKE) -C meshProcessor clean

include Makefile.deps
