FC = ifort 

CC = -traceback -p -g -check bounds -debug

all:
	$(FC) $(CC) -oHermes module_constants.f90 module_functions.f90 module_subroutines.f90 module_bonding.f90 main-iqa-prediction-stats.f90  main-microsolvation.f90  main-multipole-prediction-stats.f90  main-rdf.f90  main-structural-analysis.f90 main-Hermes.f90   

# module_bonding.f90 
# module_structure_determination.f90 

clean:
	rm *.mod
	rm *.o 
