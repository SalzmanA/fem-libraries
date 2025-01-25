#
include ../setting.mk 
#
exe=asym_elasto_damage_model asym_elasto_damage_model_AD
compile: ${exe}
asym_elasto_damage_model: asym_elasto_damage_model.cc 
	$(COMP) $(FLAGS) -o $@ $^ $(LIBS)
asym_elasto_damage_model_AD: asym_elasto_damage_model.cc 
	$(COMP) $(FLAGS) -DUSE_AD -o $@ $^ $(LIBS)
clean:
	@rm -rf ${exe}
