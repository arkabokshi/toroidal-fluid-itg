CC = mpif90
EXEC = itg
FLAGS = -O3
MAINV2 = inputdata.f90 itg.f90 deriv.f90 odesolver.f90 evolve.f90 alphainverse.f90 splines.f90

LIBS = -llapack -lblas


.PHONY: itg

itg:
	@echo "--------------------------"
	@echo "    ITG mode evolution    "
	@echo "--------------------------"
	@$(CC) -o $(EXEC) $(FLAGS) $(MAINV2) $(LIBS)
	@mpirun -np 1 nice -n 1 ./$(EXEC)

clean:
	@echo "--------------------------"
	@echo "    Cleaning up files     "
	@echo "--------------------------"
	@rm -f *.o *.mod $(EXEC) *~ *.txt *.out _tmp*



# NOTES:
# Fine as long as inputdata.f90 is being complied first!
# 'make itg' just runs the section under itg:
