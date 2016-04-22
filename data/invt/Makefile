all: wl_g wl_l nn_g nn_l

optrange="--x=0.1:0.1:50"

wl_g:
	../invtcscan.py wl_g.cfg $(optrange)

wl_l:
	../invtcscan.py wl_l.cfg $(optrange)

nn_g:
	../invtcscan.py nn_g.cfg $(optrange)

nn_l:
	../invtcscan.py nn_l.cfg $(optrange)

# quick trial run
opttrial=--x=0.1:0.2:50 --opt="--ntrials=5"

trial:
	../invtcscan.py wl_g.cfg $(opttrial)
	../invtcscan.py wl_l.cfg $(opttrial)
	../invtcscan.py nn_g.cfg $(opttrial)
	../invtcscan.py nn_l.cfg $(opttrial)
