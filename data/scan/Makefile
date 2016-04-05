all: nbsingle sigsingle sincsingle nbscan sigscan okscan

prog=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

# for the optimal alpha(t)
nbsingle: compile $(fnmdgam)
	$(prog) --nb=0.24 --samp=g --nint=10000 --fnalpha=nb0.24_g_alpha.dat
	$(prog) --nb=0.24 --samp=l --nint=10000 --fnalpha=nb0.24_l_alpha.dat
	$(prog) --nb=0.24 --samp=md --fngamma=$(fnmdgam) --nint=10000 --fnalpha=nb0.24_md_alpha.dat

sigsingle: compile $(fnmdgam)
	$(prog) --sig=10 --samp=g --nint=10000 --fnalpha=sig10_g_alpha.dat
	$(prog) --sig=10 --samp=l --nint=10000 --fnalpha=sig10_l_alpha.dat
	$(prog) --sig=10 --samp=md --fngamma=$(fnmdgam) --nint=10000 --fnalpha=sig10_md_alpha.dat


# for Figure sinc.gp
sincsingle: compile
	$(prog) --opta --okmax=8 --fnalpha=sinc_nonpbc_alpha.dat --fnwin=sinc_nonpbc_win.dat --fnwinmat=sinc_nonpbc_winmat.dat
	$(prog) --opta --okmax=4 --pbc --fnalpha=sinc_pbc_alpha.dat --fnwin=sinc_pbc_win.dat

# for the error
nbscan: compile $(fnmdgam)
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g                        --nsteps=1e7  > nbscan_t1e7_g.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l                        --nsteps=1e7  > nbscan_t1e7_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=md --fngamma=$(fnmdgam)  --nsteps=1e7  > nbscan_t1e7_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g                        --nsteps=1e8  > nbscan_t1e8_g.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l                        --nsteps=1e8  > nbscan_t1e8_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=md --fngamma=$(fnmdgam)  --nsteps=1e8  > nbscan_t1e8_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g                        --nsteps=1e9  > nbscan_t1e9_g.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l                        --nsteps=1e9  > nbscan_t1e9_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=md --fngamma=$(fnmdgam)  --nsteps=1e9  > nbscan_t1e9_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g                        --nsteps=1e10 > nbscan_t1e10_g.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l                        --nsteps=1e10 > nbscan_t1e10_l.dat
	$(prog) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=md --fngamma=$(fnmdgam)  --nsteps=1e10 > nbscan_t1e10_l.dat

sigopt=--sigdel=0.1 --sigmax=12

sigscan: compile $(fnmdgam)
	$(prog) $(sigopt) --samp=g                       --nsteps=1e8  > sigscan_t1e8_g.dat
	$(prog) $(sigopt) --samp=l                       --nsteps=1e8  > sigscan_t1e8_l.dat
	$(prog) $(sigopt) --samp=g                       --nsteps=1e10 > sigscan_t1e10_g.dat
	$(prog) $(sigopt) --samp=l                       --nsteps=1e10 > sigscan_t1e10_l.dat

okopt=--okscan=1 --okmax=99

okscan: compile $(fnmdgam)
	$(prog) $(okopt) --samp=g                       --nsteps=1e8  > okscan_t1e8_g.dat
	$(prog) $(okopt) --samp=l                       --nsteps=1e8  > okscan_t1e8_l.dat
	$(prog) $(okopt) --samp=g                       --nsteps=1e10 > okscan_t1e10_g.dat
	$(prog) $(okopt) --samp=l                       --nsteps=1e10 > okscan_t1e10_l.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

