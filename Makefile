NVCC=nvcc
NVCCFLAGS= -I. -O2 \
		   -m32 \
		   -arch sm_13 \
		   --keep --keep-dir cufiles

PTXFLAGS= -v -O2

%.o : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -c $<
