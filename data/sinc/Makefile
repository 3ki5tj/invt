all: sincsingle

prog=../../prog/predict

compile::
	make -C ../../prog predict

# for Figure sinc.gp
sincsingle: compile
	$(prog) --opta --samp=l --okmax=8       --fnwin=sinc_nonpbc_win.dat --fnwinmat=sinc_nonpbc_winmat.dat
	$(prog) --opta --samp=l --okmax=4 --pbc --fnwin=sinc_pbc_win.dat


