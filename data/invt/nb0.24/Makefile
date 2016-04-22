all: global local

global:
	../invt2run.py nb0.24_g.cfg     --c=0.5:0.5:20.0
	../invt2run.py nb0.24t1e7_g.cfg --c=0.5:0.5:20.0

local:
	../invt2run.py nb0.24_l.cfg     --c=0.5:0.5:20.0
	../invt2run.py nb0.24t1e7_l.cfg --c=0.5:0.5:20.0

# quick trial run
trial:
	../invt2run.py nb0.24_g.cfg     --c=2:2:14 --opt="--ntrials=5"
	../invt2run.py nb0.24t1e7_g.cfg --c=2:2:14 --opt="--ntrials=5"
	../invt2run.py nb0.24_l.cfg     --c=2:2:14 --opt="--ntrials=5"
	../invt2run.py nb0.24t1e7_l.cfg --c=2:2:14 --opt="--ntrials=5"
