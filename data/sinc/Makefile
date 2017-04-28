all: sincsingle

prog=../../prog/predict

compile::
	make -C ../../prog predict

# for Figure sinc.gp
sincsingle: compile
	$(prog) --okmax=10       --fnwin=sinc_nonpbc_win.dat --fnwinmat=sinc_nonpbc_winmat.dat
	$(prog) --okmax=5  --pbc --fnwin=sinc_pbc_win.dat


