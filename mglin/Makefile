main : main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o
	gcc -O3 -o main main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o
main.o : main.c nrutil.h
	gcc -c -I. main.c
mglin.o : mglin.c nrutil.h mg.h
	gcc -c -I. mglin.c
addint.o : addint.c
	gcc -c addint.c
relax.o : relax.c
	gcc -O3 -c relax.c
resid.o : resid.c
	gcc -c resid.c
rstrc.o : rstrc.c
	gcc -c rstrc.c
interp.o : interp.c
	gcc -c interp.c
slvsml.o : slvsml.c
	gcc -c slvsml.c
nrutil.o : nrutil.c
	gcc -c -I. nrutil.c
clean:
	rm -f  *.o main *~