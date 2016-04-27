all: sig

prog=../../prog/invt

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

sig: compile
	$(prog) --sig=10 --samp=g        --nequil=10000000 --ntrials=1000 > sig10_invt_g.out
	$(prog) --sig=10 --samp=g --optc --nequil=10000000 --ntrials=1000 > sig10_invtl_g.out
	$(prog) --sig=10 --samp=g --opta --nequil=10000000 --ntrials=1000 > sig10_opt_g.out
	$(prog) --sig=10 --samp=l        --nequil=10000000 --ntrials=1000 > sig10_invt_l.out
	$(prog) --sig=10 --samp=l --optc --nequil=10000000 --ntrials=1000 > sig10_invtl_l.out
	$(prog) --sig=10 --samp=l --opta --nequil=10000000 --ntrials=1000 > sig10_opt_l.out

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

