
CPP=nvcc
RM=rm -rf 
INCLUDES=-I/usr/local/cuda-7.0/include -I/home/temp/SuSparse4/include
LIBPATH64=  -L/home/temp/SuSparse4/lib -L/usr/local/cuda-7.0/lib64 -L/usr/lib64/atlas 
CFLAGS=-DGPU_BLAS -DCHOLMOD_OMP_NUM_THREADS=4   -DNPARTITION
LIBS=-lcholmod -lsuitesparseconfig -lcusolver -llapack -lcolamd -lamd -lcublas -lcusparse -lcudart -lpthread -lblas

#usage:
#	@echo "To build this example, type one of:"
#	@echo ""
#	@echo "    make build32"
#	@echo "    make build64"
#	@echo ""
#	@echo "where '32' and '64' represent the platform you wish to build for"
#	@echo ""
#	@echo "Note that this example requires Intel MKL 10.3 or greater to build"



OBJECTS= main.o\
         graph_functions.o \
	
all: $(OBJECTS) compile

compile:
	$(CPP) $(OBJECTS) -o graph $(CFLAGS) $(INCLUDES) $(LIBPATH64) $(LIBS)


main.o: main.cpp
	$(CPP) -m64 -c main.cpp $(CFLAGS) $(INCLUDES) $(LIBPATH64) $(LIBS)

graph_functions.o: graph_functions.cpp
		$(CPP) -m64 -c graph_functions.cpp $(CFLAGS) $(INCLUDES) $(LIBPATH64) $(LIBS)
	
#solver.o:
#	$(CPP) -c sources/solver.cpp 

clean:	
	$(RM) *.o	
