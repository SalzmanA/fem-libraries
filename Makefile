# compile is the default
compile:
	@cd MFEM/mechanic2d; make -f Makefile compile
	@cd FEniCSx/mechanic2d; make -f Makefile compile
clean:
	@cd MFEM/mechanic2d; make -f Makefile clean
	@cd FEniCSx/mechanic2d; make -f Makefile clean
