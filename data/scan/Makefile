all: nbscan sigscan nbsingle sigsingle

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
	$(prog) --sig=5 --samp=g --nint=10000 --fnalpha=sig5_g_alpha.dat
	$(prog) --sig=5 --samp=l --nint=10000 --fnalpha=sig5_l_alpha.dat
	$(prog) --sig=5 --samp=md --fngamma=$(fnmdgam) --nint=10000 --fnalpha=sig5_md_alpha.dat

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

sigopt=--sigdel=0.02 --sigmax=20

sigscan: compile $(fnmdgam)
	$(prog) $(sigopt) --samp=g                       --nsteps=1e8  > sigscan_t1e8_g.dat
	$(prog) $(sigopt) --samp=l                       --nsteps=1e8  > sigscan_t1e8_l.dat
	$(prog) $(sigopt) --samp=md --fngamma=$(fnmdgam) --nsteps=1e8  > sigscan_t1e8_md.dat
	$(prog) $(sigopt) --samp=g                       --nsteps=1e9  > sigscan_t1e9_g.dat
	$(prog) $(sigopt) --samp=l                       --nsteps=1e9  > sigscan_t1e9_l.dat
	$(prog) $(sigopt) --samp=md --fngamma=$(fnmdgam) --nsteps=1e9  > sigscan_t1e9_md.dat
	$(prog) $(sigopt) --samp=g                       --nsteps=1e10 > sigscan_t1e10_g.dat
	$(prog) $(sigopt) --samp=l                       --nsteps=1e10 > sigscan_t1e10_l.dat
	$(prog) $(sigopt) --samp=md --fngamma=$(fnmdgam) --nsteps=1e10 > sigscan_t1e10_md.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam): compile
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

