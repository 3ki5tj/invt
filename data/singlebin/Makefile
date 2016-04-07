all: global local

global:
	../invt2run.py singlebint1e7_g.cfg --c=0.1:0.1:5
	../invt2run.py singlebin_g.cfg     --c=0.1:0.1:5

local:
	../invt2run.py singlebint1e7_l.cfg --c=0.1:0.1:5
	../invt2run.py singlebin_l.cfg     --c=0.1:0.1:5

# quick trial run
trial:
	../invt2run.py singlebint1e7_g.cfg --c=0.5:0.5:5 --opt="--ntrials=5"
	../invt2run.py singlebin_g.cfg     --c=0.5:0.5:5 --opt="--ntrials=5"
	../invt2run.py singlebint1e7_l.cfg --c=0.5:0.5:5 --opt="--ntrials=5"
	../invt2run.py singlebin_l.cfg     --c=0.5:0.5:5 --opt="--ntrials=5"
