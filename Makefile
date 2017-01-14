IPATH = /usr/include
LPATH = /usr/lib/ -lscmpi_debug -lm
LIBS = scmpi_debug
EXES = solve tridiag_gen

all: $(EXES)

solve: solve.c cyclicreduction.c
	gcc -Wall -std=c99 -o solve solve.c cyclicreduction.c  -I$(IPATH) -L$(LPATH) -l$(LIBS)
tridiag_gen: tridiag_gen.c 
	gcc -Wall -std=c99 -o tridiag_gen tridiag_gen.c -I$(IPATH) -L$(LPATH) -l$(LIBS)

solve.o: solve.c cyclicreduction.h
	gcc -Wall -std=c99 -c solve.c cyclicreduction.h -I$(IPATH) -L$(LPATH) -l$(LIBS)
cyclicreduction.o: cyclicreduction.c cyclicreduction.h
	gcc -Wall -std=c99 -c cyclicreduction.c cyclicreduction.h -I$(IPATH) -L$(LPATH) -l$(LIBS)

clean:
	\rm -f *.o $(EXES) PI* *gch tridiag*txt

