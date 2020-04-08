FC=gfortran
F2PY=f2py
HELM_TABLE_DIR=$(shell pwd)
HELM_TABLE_NAME=helm_table.dat

all: module

module: helmholtz.o eosfxt.o
	${F2PY} -m fhelmholtz --fcompiler=${FC} -c pycall.f90 -I helmholtz.o
	${F2PY} -m ftimmes --fcompiler=${FC} -c pycall_eosfxt.f90 -I eosfxt.o

test: test.o helmholtz.o
	${FC} -o test.x test.o helmholtz.o
	./test.x

helmholtz.o: helmholtz.f90 const.dek implno.dek vector_eos.dek
	${FC} -cpp -DTBLPATH="'${HELM_TABLE_DIR}/${HELM_TABLE_NAME}'" -ffree-line-length-none -c -fPIC $<

eosfxt.o: eosfxt.f90 const.dek implno.dek vector_eos.dek
	${FC} -c -fPIC $<

%.o : %.f90
	${FC} -c $<

clean:
	rm -f *.o *.so *.x
