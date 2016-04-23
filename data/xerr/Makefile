all: sig

prog=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

sig: compile $(fnmdgam)
	$(prog) --sig=10 --samp=g                        --nint=100000 --fnxerr=xerr_g.dat
	$(prog) --sig=10 --samp=l                        --nint=100000 --fnxerr=xerr_l.dat
	$(prog) --sig=10 --samp=md --fngamma=$(fnmdgam)  --nint=100000 --fnxerr=xerr_md.dat

sinc: compile $(fnmdgam)
	$(prog) --okmax=4 --samp=g                       --nint=100000 --fnxerr=xerr_g.dat
	$(prog) --okmax=4 --samp=l                       --nint=100000 --fnxerr=xerr_l.dat
	$(prog) --okmax=4 --samp=md --fngamma=$(fnmdgam) --nint=100000 --fnxerr=xerr_md.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

