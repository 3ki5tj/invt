all: gamma

pred=../../prog/predict

fn1stepgam=1stepgamma.dat
fnmdgam=mdgamma.dat
fnlgam=lgamma.dat

compile::
	make -C ../../prog

gamma: compile $(fn1stepgam) $(fnmdgam) $(fnlgam)
	paste $(fn1stepgam) $(fnmdgam) $(fnlgam) | tail -n +3 > gamma.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$@ --gamnsteps=10000000000 --try=0

# use `invt` to generate lgamma.dat
$(fnlgam):
	../../prog/invt --samp=l  --gam --fngamma=$@ --gamnsteps=10000000000 --try=0
	
$(fn1stepgam):
	$(pred) --samp=l --fngamma=$(fn1stepgam)
