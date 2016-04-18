all: gamma

pred=../../prog/predict

fnlgam=lgamma.dat
fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

gamma: compile $(fnmdgam)
	$(pred) --samp=l --fngamma=$(fnlgam)
	paste $(fnlgam) $(fnmdgam) | tail -n +3 > gamma.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

