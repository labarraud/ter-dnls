FC = gfortran
FCFLAGS=-g -fcheck=all -Warray-bounds

PROG=ter
SOURCES=  resolveLUComplex.o creesys.o exact.o fonction_u.o equation.o resolve.o tableau.o ecrirefichier.o erreur.o main.o

all:$(PROG)

$(PROG):$(SOURCES)
	$(FC) $(FCFLAGS) -o $@ $^

%.o:%.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod

cleanall:
	rm -f $(PROG) *.o *.mod
