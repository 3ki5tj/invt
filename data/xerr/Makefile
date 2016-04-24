all: sig

prog=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

sig: compile $(fnmdgam)
	$(prog) --sig=10 --nsteps=1e8   --samp=g                        --nint=100000 --fnxerr=xerr_t1e8_g.dat
	$(prog) --sig=10 --nsteps=1e8   --samp=l                        --nint=100000 --fnxerr=xerr_t1e8_l.dat
	$(prog) --sig=10 --nsteps=1e8   --samp=md --fngamma=$(fnmdgam)  --nint=100000 --fnxerr=xerr_t1e8_md.dat
	$(prog) --sig=10 --nsteps=1e9   --samp=g                        --nint=100000 --fnxerr=xerr_t1e9_g.dat
	$(prog) --sig=10 --nsteps=1e9   --samp=l                        --nint=100000 --fnxerr=xerr_t1e9_l.dat
	$(prog) --sig=10 --nsteps=1e9   --samp=md --fngamma=$(fnmdgam)  --nint=100000 --fnxerr=xerr_t1e9_md.dat
	$(prog) --sig=10 --nsteps=1e10  --samp=g                        --nint=100000 --fnxerr=xerr_t1e10_g.dat
	$(prog) --sig=10 --nsteps=1e10  --samp=l                        --nint=100000 --fnxerr=xerr_t1e10_l.dat
	$(prog) --sig=10 --nsteps=1e10  --samp=md --fngamma=$(fnmdgam)  --nint=100000 --fnxerr=xerr_t1e10_md.dat
	$(prog) --sig=10 --nsteps=1e11  --samp=g                        --nint=100000 --fnxerr=xerr_t1e11_g.dat
	$(prog) --sig=10 --nsteps=1e11  --samp=l                        --nint=100000 --fnxerr=xerr_t1e11_l.dat
	$(prog) --sig=10 --nsteps=1e11  --samp=md --fngamma=$(fnmdgam)  --nint=100000 --fnxerr=xerr_t1e11_md.dat

sinc: compile $(fnmdgam)
	$(prog) --okmax=4 --samp=g                       --nint=100000 --fnxerr=xerr_g.dat
	$(prog) --okmax=4 --samp=l                       --nint=100000 --fnxerr=xerr_l.dat
	$(prog) --okmax=4 --samp=md --fngamma=$(fnmdgam) --nint=100000 --fnxerr=xerr_md.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

