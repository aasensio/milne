# Compiler 
# ifort -> Intel Fortran
# gfortran -> gfortran

# The following lines are used to choose the compiler
# We have tested gfortran and ifort
# In case you want to use another compiler, use equivalent
# keywords
#  FLAGS -> is used to indicate that the preprocessor has to be invoked
#  OPTIONS -> these are general compilation flags that, for the moment, only indicate
#             that the object file is generated, without linking
# It should be easy to find the equivalent flags in your compiler

# gfortran (comment these lines if you want to use another compiler)
# COMPILER = gfortran
# PYCOMPILER = gnu95
# FLAGS =
# OPTIONS = -c -ffree-line-length-none -fPIC
 
# ifort (comment these lines if you want to use another compiler)
COMPILER = ifort
PYCOMPILER = intel
LIBS = -lifcore -lsvml
OPTIONS = -c -fPIC -O3
 
FINAL_EXECUTABLE = milne
F2PY = /usr/pkg/python/Python-2.7.4/bin/f2py

# 

all:
	make pymilne
		
pymilne: vars.o maths.o atomic.o milne.o
	$(F2PY) $(LIBS) -c milne.pyf vars.o maths.o atomic.o milne.o

clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "milne.so" -delete ; find . -maxdepth 1 -name "*.pyf" -delete
	find ../ -maxdepth 2 -name "milne.so" -delete
	
vars.o: vars.f90
	$(COMPILER) $(OPTIONS)  vars.f90

maths.o: maths.f90 vars.o
	$(COMPILER) $(OPTIONS)  maths.f90

atomic.o: atomic.f90 maths.o vars.o
	$(COMPILER) $(OPTIONS)  atomic.f90
	
milne.o: milne.f90
	$(COMPILER) $(OPTIONS)  milne.f90
	$(F2PY) --overwrite-signature -m milne -h milne.pyf milne.f90
