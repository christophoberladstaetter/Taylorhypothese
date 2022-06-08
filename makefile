exec = vor2
objs = vor2.o 

#LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lpthread -lfftw3_threads -lfftw3 #-pg
####LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lfftw3_omp -lfftw3
LFLAGS = -mcmodel=large -march=native -lm -O2 -fopenmp -lfftw3_omp -lfftw3 

CFLAGS = -I.
LIBS = /usr/lib64 

all: $(exec)

$(exec): $(objs)
	g++  $(exec).C -o $@  $(LFLAGS) -L$(LIBS)
	rm *.o

ghw.o: vor2.C 

	g++  -c $< -o $@

.PHONY : clean

clean :
	rm -f $(exec) $(objs) *~ 
