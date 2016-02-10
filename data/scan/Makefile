all: nbscan sigscan nbsingle sigsingle

prog=../../prog/predict

compile::
	make -C ../../prog predict

# for the optimal alpha(t)
nbsingle: compile
	$(prog) --nb=0.24 --samp=g --fnalpha=nb0.24_g_alpha.dat
	$(prog) --nb=0.24 --samp=l --fnalpha=nb0.24_l_alpha.dat

sigsingle: compile
	$(prog) --sig=5 --samp=g --fnalpha=sig5_g_alpha.dat
	$(prog) --sig=5 --samp=l --fnalpha=sig5_l_alpha.dat

# for the error
nbscan: compile
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g > nbscan_g.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l > nbscan_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g --nsteps=1e7 > nbscan_t1e7_g.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l --nsteps=1e7 > nbscan_t1e7_l.dat

sigscan: compile
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e6  > sigscan_t1e6_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e6  > sigscan_t1e6_l.dat
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e7  > sigscan_t1e7_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e7  > sigscan_t1e7_l.dat
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e8  > sigscan_t1e8_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e8  > sigscan_t1e8_l.dat
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e9  > sigscan_t1e9_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e9  > sigscan_t1e9_l.dat
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e10 > sigscan_t1e10_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e10 > sigscan_t1e10_l.dat
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e11 > sigscan_t1e11_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e11 > sigscan_t1e11_l.dat
	$(prog) --sigdel=0.02 --samp=g --nsteps=1e12 > sigscan_t1e12_g.dat
	$(prog) --sigdel=0.02 --samp=l --nsteps=1e12 > sigscan_t1e12_l.dat

