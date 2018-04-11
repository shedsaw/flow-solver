#include folder of headers
#INCLUDE_PATH = -I/home/ssawyer/RESEARCH/FIELD_SOLVER_2D/include -I/usr/local/gsl/include
INCLUDE_PATH = -I/home/ssawyer/RESEARCH/FIELD_SOLVER_2D/include -I/home/ssawyer/lib/include/

#LIB_PATH = -L/usr/local/gsl/lib
LIB_PATH = -L/home/ssawyer/lib/lib

# Standard compiler
CC = gcc

# optimized and correct code for gcc (pedantic not supported in icc)
# CFLAGS = -Wno-deprecated -Wall -ansi -pedantic -O3

# icc
#  CFLAGS = -O3

# code needing debugging
#  CFLAGS = -g -Wno-deprecated -Wall -ansi
  CFLAGS = -Wno-deprecated -Wall -ansi -pedantic -O3

#math library, std c++ lib
  LIBS = -lm -lstdc++ -lblas -llapack -lgsl -lgslcblas
#  LIBS = -lm -lstdc++ -lblas -llapack -lgsl -lgslcblas -lefence 

OBJ_OUT = main.o cv_calc.o geometry.o grid_io.o linear_algebra.o maps.o \
	mesh_io.o mesh_utils.o metrics.o params.o reorder.o gradient.o \
	hessian.o flux.o residual.o init_memory.o ic.o limiter.o solve.o cvbc.o \
	jacobian.o reconstruct.o curved_boundaries.o compute_derivatives.o

all: field_solver

field_solver: $(OBJ_OUT)
	$(CC) $(C_FLAGS) $(LIB_PATH) -o $@ $(OBJ_OUT) $(LIBS)

clean:
	\rm -rf *.o field_solver

#.C.o:
#	$(CC) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) -c $<

.C.o:
	$(CC) $(CFLAGS) $(INCLUDE_PATH) $(LIB_PATH) $(LIBS) -c $<
