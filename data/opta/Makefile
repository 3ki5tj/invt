all: wlsingle nbsingle sigsingle

prog=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

# for the optimal alpha(t)
wlsingle: compile $(fnmdgam)
	$(prog) --samp=g                       --nint=10000 --aresamp --fnalpha=wl_g_alpha.dat
	$(prog) --samp=l                       --nint=10000 --aresamp --fnalpha=wl_l_alpha.dat
	$(prog) --samp=md --fngamma=$(fnmdgam) --nint=10000 --aresamp --fnalpha=wl_md_alpha.dat

nbsingle: compile $(fnmdgam)
	$(prog) --nb=0.24 --samp=g                       --nint=10000 --aresamp --fnalpha=nb0.24_g_alpha.dat
	$(prog) --nb=0.24 --samp=l                       --nint=10000 --aresamp --fnalpha=nb0.24_l_alpha.dat
	$(prog) --nb=0.24 --samp=md --fngamma=$(fnmdgam) --nint=10000 --aresamp --fnalpha=nb0.24_md_alpha.dat

sigsingle: compile $(fnmdgam)
	$(prog) --sig=10 --samp=g                       --nint=10000 --aresamp                             --fnalpha=sig10_g_alpha.dat
	$(prog) --sig=10 --samp=l                       --nint=10000 --aresamp                             --fnalpha=sig10_l_alpha.dat
	$(prog) --sig=10 --samp=md --fngamma=$(fnmdgam) --nint=10000 --aresamp                             --fnalpha=sig10_md_alpha.dat
	$(prog) --sig=10 --samp=g                       --nint=10000 --aresamp --nsteps=50000000 --a0=2e-4 --fnalpha=sig10_g_t5e7_alpha.dat
	$(prog) --sig=10 --samp=l                       --nint=10000 --aresamp --nsteps=50000000 --a0=2e-4 --fnalpha=sig10_l_t5e7_alpha.dat
	$(prog) --sig=10 --samp=md --fngamma=$(fnmdgam) --nint=10000 --aresamp --nsteps=50000000 --a0=2e-4 --fnalpha=sig10_md_t5e7_alpha.dat
	$(prog) --sig=10 --samp=g                       --nint=10000                                       --fnalpha=sig10_g_alpha_q.dat
	$(prog) --sig=10 --samp=l                       --nint=10000                                       --fnalpha=sig10_l_alpha_q.dat
	$(prog) --sig=10 --samp=md --fngamma=$(fnmdgam) --nint=10000                                       --fnalpha=sig10_md_alpha_q.dat
	$(prog) --sig=10 --samp=g                       --nint=10000           --nsteps=50000000 --a0=2e-4 --fnalpha=sig10_g_t5e7_alpha_q.dat
	$(prog) --sig=10 --samp=l                       --nint=10000           --nsteps=50000000 --a0=2e-4 --fnalpha=sig10_l_t5e7_alpha_q.dat
	$(prog) --sig=10 --samp=md --fngamma=$(fnmdgam) --nint=10000           --nsteps=50000000 --a0=2e-4 --fnalpha=sig10_md_t5e7_alpha_q.dat


# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

