CC=gcc
MPICC=mpicc
CFLAGS=-Wall -g -std=gnu99 -O3

.PHONY: exe clean realclean


# === Executables

exe: path-omp.x path-mpi.x

path-omp.x: path-omp.o mt19937p.o
	$(CC) -fopenmp $(CFLAGS) $^ -o $@

path-omp.o: path-omp.c
	$(CC) -c -fopenmp $(CFLAGS) $<

path-mpi.x: path-mpi.o mt19937p.o
	$(MPICC) $(CFLAGS) $^ -o $@

path-mpi.o: path-mpi.c
	$(MPICC) -c $(CFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<

# == Plots
OMPRESULTS = exp_omp_1 exp_omp_2 exp_omp_3 exp_omp_4 \
exp_omp_5 exp_omp_6 exp_omp_7 exp_omp_8

omp_plots:
	python omp_plot.py strong $(OMPRESULTS) omp_strong.pdf
	python omp_plot.py weak $(OMPRESULTS) omp_weak.pdf

# === Documentation

main.pdf: main.tex path-omp.tex

path-omp.tex: path-omp.c
	dsbweb -o $@ -c $^

%.pdf: %.tex
	pdflatex $<
	pdflatex $<


# === Cleanup and tarball

clean:
	rm -f *.o
	rm -f main.aux main.log main.out

realclean: clean
	rm -f path-omp.x path-mpi.x path-omp.tex main.pdf *.*~ run_omp*.qsub* *.pyc

tgz:
	( cd ..; tar -czf path.tgz \
	    path/Makefile path/main.tex \
	    path/path-omp.c \
	    path/run-omp.qsub path/run-mpi.qsub \
	    path/mt19937p.h path/mt19937p.c )
