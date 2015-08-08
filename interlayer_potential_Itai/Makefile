#linux_LIB = -lblas  -llapack -lg2c -lblas -lfftw
linux_LIB = -L/usr/lib64/ 
linux_MYLIB = -L${HOME}/lib/
linux_CLIB =  -lm -lgcc -lstdc++
linux_FCLIB = 	
OSF1_LIB = -ldxml  /usr/local/lib/lapack_alpha.a
OSF1_MYLIB =
OSF1_CLIB =  -lm -lfor
OSF1_FCLIB =
IRIX_LIB =  -lcomplib.sgimath -L/usr/local/lib64/nag17 -lnag_n64_m4
IRIX_MYLIB = 
IRIX_CLIB =  -lm
IRIX_FCLIB = 	

LIB =  $(${OSTYPE}_MYLIB) $(${OSTYPE}_LIB) $(${OSTYPE}_CLIB) $(${OSTYPE}_FCLIB) $(${OSTYPE}_MYLIB)

# flags ...
linux_OPTFLAGS = -O4 -fopenmp -ffast-math -Wno-deprecated
linux_CFLAGS = -DTCPUNIX $(${OSTYPE}_OPTFLAGS)
linux_FFLAGS = $(${OSTYPE}_OPTFLAGS)

OSF1_OPTFLAGS = -O5
OSF1_CFLAGS = $(${OSTYPE}_OPTFLAGS) 
OSF1_FFLAGS = -no_main $(${OSTYPE}_OPTFLAGS)

IRIX_OPTFLAGS = -64 -O4 -mips4 -OPT:Olimit=0
IRIX_CFLAGS = $(${OSTYPE}_OPTFLAGS) 
IRIX_FFLAGS = $(${OSTYPE}_OPTFLAGS) 

#LF_CFLAGS = -fopenmp
#LF_FFLAGS = -fopenmp

DIAG = diag
QRFAC = qrfac
MAINNAM = 27.5.14_MD

#compiler ...
linux_CC = g++
linux_FF = g77
linux_LD = g++

IRIX_CC = CC
IRIX_FF = f77
IRIX_LD = CC

OSF1_CC = c++
OSF1_FF = f77
OSF1_LD = c++

OBJECTS = \
main.o ReadAtoms.o Bond_Parameters.o read_parameters.o Potential.o calc_rij.o Propagate.o Print.o rand.o  Vel.o Print_to_log.o Quench.o gauss-rand.o
# compilation ...

# (use gmake for IRIX).
.f.o:
	$(${OSTYPE}_FF) $(${OSTYPE}_FFLAGS) -c  $*.f
.c.o:
	$(${OSTYPE}_CC) $(${OSTYPE}_CFLAGS) -c  $*.c  
 
$(MAINNAM): $(OBJECTS) 
	$(${OSTYPE}_LD) -o $(MAINNAM).x $(${OSTYPE}_CFLAGS) $(OBJECTS) $(LIB)

$(DIAG).o:	$(DIAG).f
$(QRFAC).o:	$(QRFAC).f

clean:
	/bin/rm *.o *.x
