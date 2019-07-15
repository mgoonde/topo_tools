F90 = gfortran
CC  = gcc

FFLAGS = -Wall -O0 -fdefault-real-8 -cpp -fcheck=all

NAUT = ./nauty26r11/




default: coords_hash.x




ffnauty.o:
	$(CC) -c ffnauty.c -I./$(NAUT)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

f90nautyinterf.o: f90nautyinterf.f90
	$(F90) $(FFLAGS) -c $<








coords_hash.x: f90nautyinterf.o basic_module.o graph_module.o coords_hash.o ffnauty.o
	$(F90) $(FFLAGS) -o coords_hash.x $^ $(NAUT)/nauty.a



clean:
	rm -rf *.o *.mod *.x                                     
