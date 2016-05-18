iarundat = iarun_sig10_g.dat iarun_sig10_l.dat

all: $(iarundat)

$(iarundat): % : Makefile
	../scanave.py ../iascan[1-9]*/$@ -o $@
