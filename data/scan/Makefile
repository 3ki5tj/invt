all: sigprd okprd sigrun okrun

pred=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

# for the error
nbprd: compile $(fnmdgam)
	$(pred) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=g                        --nsteps=1e8  > nbscan_t1e8_g.dat
	$(pred) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=l                        --nsteps=1e8  > nbscan_t1e8_l.dat
	$(pred) --nbmin=-1.0 --nbdel=0.001 --nbmax=0.25 --samp=md --fngamma=$(fnmdgam)  --nsteps=1e8  > nbscan_t1e8_l.dat

sigprdopt=--sigdel=0.1 --sigmax=12

sigprd: compile $(fnmdgam)
	$(pred) $(sigprdopt) --samp=g                       --nsteps=1e8  > sigprd_t1e8_g.dat
	$(pred) $(sigprdopt) --samp=l                       --nsteps=1e8  > sigprd_t1e8_l.dat
	$(pred) $(sigprdopt) --samp=g                       --nsteps=1e11 > sigprd_t1e11_g.dat
	$(pred) $(sigprdopt) --samp=l                       --nsteps=1e11 > sigprd_t1e11_l.dat

sigrunopt=--sig=0.5:0.5:12

sigrun: compile $(fnmdgam)
	../invtscan.py $(sigrunopt) global.cfg -o sigrun_t1e8_g.dat
	../invtscan.py $(sigrunopt) local.cfg  -o sigrun_t1e8_l.dat

okprdopt=--okscan=1 --okmax=99

okprd: compile $(fnmdgam)
	$(pred) $(okprdopt) --samp=g                       --nsteps=1e8  > okprd_t1e8_g.dat
	$(pred) $(okprdopt) --samp=l                       --nsteps=1e8  > okprd_t1e8_l.dat
	$(pred) $(okprdopt) --samp=g                       --nsteps=1e11 > okprd_t1e11_g.dat
	$(pred) $(okprdopt) --samp=l                       --nsteps=1e11 > okprd_t1e11_l.dat

okrunopt=--ok=1:1:20

okrun: compile $(fnmdgam)
	../invtscan.py $(okrunopt) global.cfg -o okrun_t1e8_g.dat
	../invtscan.py $(okrunopt) local.cfg  -o okrun_t1e8_l.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

