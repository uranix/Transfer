CUDA_INSTALL_PATH=/usr/local/cuda
NVCC=$(CUDA_INSTALL_PATH)/bin/nvcc
NVCCFLAGS= -I. -O2 \
		   -m32 \
		   -arch sm_13 \
		   --keep --keep-dir cufiles
CFLAGS=-m32 -O2
CXXFLAGS=-m32 -O2
LDFLAGS=-m32
PTXFLAGS= -v -O2
LIBS= -L/usr/local/cuda/lib/ -lcudart

.PHONY: all clean

all: main

%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -c $<

main: main.o wrapper.o 
	gcc $(LDFLAGS) -o $@ $^ $(LIBS)

clean: 
	rm -f *.o
