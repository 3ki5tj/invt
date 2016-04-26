all: iaprd iarun_l iarun_g

pred=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

iaprdopt=--iamin=5e-7 --iamax=2e-2 --iadel=0.01

iaprd: compile $(fnmdgam)
	$(pred) $(iaprdopt) --sig=10 --samp=l                       > iaprd_sig10_l.dat
	$(pred) $(iaprdopt) --sig=10 --samp=g                       > iaprd_sig10_g.dat
	$(pred) $(iaprdopt) --sig=10 --samp=md --fngamma=$(fnmdgam) > iaprd_sig10_md.dat
	$(pred) $(iaprdopt) --sig=3  --samp=l                       > iaprd_sig3_l.dat
	$(pred) $(iaprdopt) --sig=3  --samp=g                       > iaprd_sig3_g.dat
	$(pred) $(iaprdopt) --sig=3  --samp=md --fngamma=$(fnmdgam) > iaprd_sig3_md.dat

iarunopt=--ia=5e-7:0.2:2e-2

iarun_l: compile
	../invtscan.py $(iarunopt) --opt="--sig=10" local.cfg  -o iarun_sig10_l.dat
	../invtscan.py $(iarunopt) --opt="--sig=3"  local.cfg  -o iarun_sig3_l.dat

iarun_g: compile
	../invtscan.py $(iarunopt) --opt="--sig=10" global.cfg  -o iarun_sig10_g.dat
	../invtscan.py $(iarunopt) --opt="--sig=3"  glocal.cfg  -o iarun_sig3_g.dat

iarun_trial: compile
	../invtscan.py $(iarunopt) --opt="--sig=10 --ntrials=3" local.cfg   -o iarun_sig10_l.dat
	../invtscan.py $(iarunopt) --opt="--sig=10 --ntrials=3" global.cfg  -o iarun_sig10_g.dat

# use `invt` to generate mdgamma.dat
$(fnmdgam):
	../../prog/invt --samp=md --gam --fngamma=$(fnmdgam) --gamnsteps=10000000000 --try=0

