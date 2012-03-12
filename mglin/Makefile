main : main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o 
	gcc -O3 -o main main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o
main.o : main.c nrutil.h
	gcc -c -g -I. main.c
mglin.o : mglin.c nrutil.h mg.h
	gcc -c -g -I. mglin.c
addint.o : addint.c
	gcc -c -g addint.c
relax.o : relax.c
	gcc -O3 -c -g relax.c
resid.o : resid.c
	gcc -c -g resid.c
rstrc.o : rstrc.c
	gcc -c -g rstrc.c
interp.o : interp.c
	gcc -c -g interp.c
slvsml.o : slvsml.c
	gcc -c -g slvsml.c
nrutil.o : nrutil.c
	gcc -c -g nrutil.c
clean:
	rm -f  *.o main *~
