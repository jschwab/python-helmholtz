FC=gfortran
F2PY=f2py
HELM_TABLE_DIR=$(shell pwd)
HELM_TABLE_NAME=helm_table.dat

all: module

module: helmholtz.o
	${F2PY} -m fhelmholtz --fcompiler=${FC} -c pycall.f90 -I helmholtz.o

test: test.o helmholtz.o
	${FC} -o test.x test.o helmholtz.o
	./test.x

helmholtz.o : helmholtz.f90
	${FC} -cpp -DTBLPATH="'${HELM_TABLE_DIR}/${HELM_TABLE_NAME}'" -ffree-line-length-none -c -fPIC $<

%.o : %.f90
	${FC} -c $<

clean:
	rm -f *.o *.x
