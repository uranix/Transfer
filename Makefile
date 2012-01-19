NVCC=nvcc
NVCCFLAGS= -I. -O2 \
		   -ccbin "D:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin" \
		   -arch sm_13 \
		   --keep --keep-dir cufiles

PTXFLAGS= -v -O2

%.obj : %.cu
	$(NVCC) $(NVCCFLAGS) -Xptxas "$(PTXFLAGS)" -c $<
