FC=gfortran

module: helmholtz.o
	f2py -m fhelmholtz --fcompiler=${FC} -c pycall.f90 -I helmholtz.o

test: test.o helmholtz.o
	${FC} -o test.x test.o helmholtz.o
	./test.x

%.o : %.f90
	${FC} -c -fPIC $<

clean:
	rm -f *.o *.x
