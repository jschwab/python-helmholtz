FC=gfortran

test: test.o helmholtz.o
	${FC} -o test.x test.o helmholtz.o
	./test.x

%.o : %.f90
	${FC} -c $<

clean:
	rm -f *.o *.x
