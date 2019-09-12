###################################
# Build and install MT3DNC        #
# By Zhuxiaoxiong                 #
###################################

CC   =/vol7/home/test653/zxx/mpi3-gcc/bin/mpicc
CXX  =/vol7/home/test653/zxx/mpi3-gcc/bin/mpicxx


LIB_SuperLU =/vol7/home/test653/zxx/lib/libsuperlu_dist.a
LIB_OpenBLAS=/vol7/home/test653/zxx/lib/libopenblas.a
LIB_parmetis=/vol7/home/test653/zxx/lib/libparmetis.a
LIB_metis   =/vol7/home/test653/zxx/lib/libmetis.a
LIB_mpi     =/vol7/home/test653/zxx/mpi3-gcc/lib/libmpi.a \
             /vol7/home/test653/zxx/mpi3-gcc/lib/libmpicxx.a

INCLUDE =-I./include
CFLAGES = -c -g -O2 -fopenmp 
LDFLAGES= -g -O2 -fopenmp -o

target=MTSZ
#obj   =main.o prepocess.o global_varible.o post_prepocess.o \
       FEM_compute.o para.o solver.o
obj  = *.o
#src   =main.c prepocess.c global_varible.c post_prepocess.c \
       FEM_compute.c para.c solver.c
src   = *.c
$(target):$(obj) $(LIB_SuperLU) 
#$(LIB_OpenBLAS) $(LIB_parmetis)$(LIB_metis) $(LIB_mpi)
	$(CXX) $(LDFLAGES) $(target) $(obj) $(LIB_SuperLU) $(LIB_OpenBLAS) \
	$(LIB_parmetis) $(LIB_metis) -lm

$(obj):$(src)
	$(CC) $(CFLAGES) $(src) $(INCLUDE)

clean: 
	rm $(obj) $(target)
