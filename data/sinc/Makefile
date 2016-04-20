all: sincsingle

prog=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

# for Figure sinc.gp
sincsingle: compile
	$(prog) --opta --samp=l --okmax=6       --aresamp --fnalpha=sinc_nonpbc_alpha.dat --fnwin=sinc_nonpbc_win.dat --fnwinmat=sinc_nonpbc_winmat.dat
	$(prog) --opta --samp=l --okmax=3 --pbc --aresamp --fnalpha=sinc_pbc_alpha.dat --fnwin=sinc_pbc_win.dat


