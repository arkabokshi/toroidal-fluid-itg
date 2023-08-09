CC = gfortran #mpif90
EXEC = itg
#FLAGS = -O3 - -g -Og -fcheck=all 
FLAGS =-fdefault-real-8 -fdefault-double-8 -fimplicit-none -O3 -fopenmp # -fcheck=all,no-array-temps -ffpe-trap=invalid,zero -Og -g -Wall -fimplicit-none #-fsanitize=address,undefined -fno-omit-frame-pointer -fimplicit-none
MAINV2 = inputdata.f90 lapack_wrap.f90 deriv.f90 alphainverse.f90 evolve.f90 odesolver.f90 itg.f90

LIBS = -llapack -lblas


.PHONY: itg

itg:
	@echo "--------------------------"
	@echo "    ITG mode evolution    "
	@echo "--------------------------"
	@$(CC) -o $(EXEC) $(FLAGS) $(MAINV2) $(LIBS)
	@./$(EXEC)

clean:
	@echo "--------------------------"
	@echo "    Cleaning up files     "
	@echo "--------------------------"
	@rm -f *.o *.mod $(EXEC) *~ *.txt *.out _tmp*



# NOTES:
# Fine as long as inputdata.f90 is being complied first!
# 'make itg' just runs the section under itg:
