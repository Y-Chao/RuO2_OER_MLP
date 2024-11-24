source_1 = arrayprocess.f90 structure.f90 bondmatrix.f90
source_2 = fakebondmatrix.f90 bond_func.f90 

#FC = ifort
FC = gfortran
FCFLAGS = -shared -fPIC -O2 -std=f2008 
#FCFLAGS = -shared -fPIC -O2 -xHost -check bounds

# Default target: build both shared libraries
all: bondmatrix.so checkminbond.so

bondmatrix.so: ${source_1}
	${FC} ${FCFLAGS} -o bondmatrix.so ${source_1}

checkminbond.so: ${source_2}
	${FC} ${FCFLAGS} -o checkminbond.so ${source_2}

.PHONY: clean

clean:
	rm -f *.o *.so *.mod
