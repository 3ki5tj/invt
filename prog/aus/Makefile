CC = icc -Wall -Wremarks -O2 -g
LM =
#CC = gcc -O3 -Wall -Wextra -g
#LM = -lm

CCd = gcc -Wall -Wextra -O0 -g
LMd = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h ../*.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) \
	  a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos *.log *.his vgcore*
	rstrip.py -Rlv

Bossman: clean
	rsync -avz --exclude=bak * ../data/*.py /Bossman/cz1/invt/prog/

Dropbox: clean
	$(MAKE) -C .. $@

