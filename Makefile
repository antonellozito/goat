# A very simple and straightforward makefile without too many options 
# yet. To be improved in the future. 

# Variables
FC = gfortran # compiler 
CFLAGS = -c -g -Og -Wall # compiler flags

# Linking
a.out: gdmod_types.o gdmod_userinput.o MainRunFileGridDeformation.o
	$(FC) MainRunFileGridDeformation.o gdmod_types.o gdmod_userinput.o

#MainRunFileGridDeformation.o: src/Modules/gdmod_types.o

# Compiling
gdmod_types.o: src/Modules/gdmod_types.F90
	$(FC) $(CFLAGS) src/Modules/gdmod_types.F90

gdmod_userinput.o: src/Modules/gdmod_userinput.F90
	$(FC) $(CFLAGS) src/Modules/gdmod_userinput.F90

MainRunFileGridDeformation.o: Runfiles/MainRunFileGridDeformation.F90
	$(FC) $(CFLAGS) Runfiles/MainRunFileGridDeformation.F90


# Cleanup
.PHONY: clean
clean: 
	rm *.o a.out

# Run
.PHONY: gd
gd: 
	a.out
	make clean