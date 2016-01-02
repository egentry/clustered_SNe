
.PHONY: clean

spherical_standard_omega0p00_nstep_01000.dat : sedov3
	./sedov3 1000

sedov3 : sedov3.f
	gfortran sedov3.f -o sedov3

clean:
	rm -f sedov3 *.dat *.o