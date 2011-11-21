FLAGS = -O3 -threaded -rtsopts
PROFILE = -prof -auto-all -eventlog
HPC = -fhpc
FILES = repl.hs eval.hs parse.hs types.hs dispatch.hs keywords.hs operators.hs

all :
	ghc $(FLAGS) -o scheme522 $(FILES)

profile :
	ghc $(FLAGS) $(PROFILE) -o scheme522 $(FILES)

hpc :
	ghc $(FLAGS) $(HPC) -o scheme522 $(FILES)

clean :
	rm -f *.o *~ *.hi scheme522
