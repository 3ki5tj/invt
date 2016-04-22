all: sincsingle

prog=../../prog/predict

fnmdgam=mdgamma.dat

compile::
	make -C ../../prog

# for Figure sinc.gp
sincsingle: compile
	$(prog) --opta --samp=l --okmax=4       --fnwin=sinc_nonpbc_win.dat --fnwinmat=sinc_nonpbc_winmat.dat
	$(prog) --opta --samp=l --okmax=2 --pbc --fnwin=sinc_pbc_win.dat


