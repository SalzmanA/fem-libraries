#
include ../setting.mk 
# core rules
exe= asym_elasto_damage_model_symb_sym asym_elasto_damage_model_symb  asym_elasto_damage_model_manual asym_elasto_damage_model_ufl
compile: ${exe}
asym_elasto_damage_model_symb_sym: asym_elasto_damage_model.cc asym_symb_sym.o
	$(COMP) $(FLAGS) -DSYMB_SYM -o $@ $^ $(LIBS)
asym_elasto_damage_model_symb: asym_elasto_damage_model.cc asym_symb.o
	$(COMP) $(FLAGS) -DSYMB -o $@ $^ $(LIBS)
asym_elasto_damage_model_manual: asym_elasto_damage_model.cc asym_manual.o
	$(COMP) $(FLAGS) -DMANUAL -o $@ $^ $(LIBS)
asym_elasto_damage_model_ufl: asym_elasto_damage_model.cc asym_ufl.o
	$(COMP) $(FLAGS) -DUFL_POTENTIAL -o $@ $^ $(LIBS)
%.o:%.c
	gcc -c ${OPT} -o $@ $< -I${UFCX_DIR}
%.c:%.py
	ffcx $<
.SECONDARY: asym_symb_sym.c asym_symb.c asym_manual.c asym_ufl.c
