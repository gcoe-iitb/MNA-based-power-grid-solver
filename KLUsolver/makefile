
CUDAROOT=/usr/local/cuda-6.5
SUITESPARSEROOT=/usr/local
CPP=g++
INCLUDES=-I$(CUDAROOT)/include -I/usr/include
LIBPATH64=-L$(SUITESPARSEROOT)/lib -L$(CUDAROOT)/lib64
CFLAGS= -O3 -DGPU_BLAS -DCHOLMOD_OMP_NUM_THREADS=4 -DNPARTITION -DNPARTITION -fopenmp -g
LIBS= -lklu -lbtf -lcholmod -lsuitesparseconfig -llapack -lcolamd -lamd -lcublas -lcudart -lcusparse -lpthread -lblas

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


main.o: sources/main.cpp
	$(CPP) -m64 -c sources/main.cpp $(CFLAGS) $(INCLUDES) $(LIBPATH64) $(LIBS)

graph_functions.o: sources/graph_functions.cpp
		$(CPP) -m64 -c sources/graph_functions.cpp $(CFLAGS) $(INCLUDES) $(LIBPATH64) $(LIBS)
	
#solver.o:
#	$(CPP) -c sources/solver.cpp 

clean:	
	rm *.o	
