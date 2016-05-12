CC = icc -Wall -Wremarks -O2 -g
LM =
#CC = gcc -O3 -Wall -Wextra -g
#LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) \
	  a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos *.log
	rstrip.py -Rlv

Bossman: all
	rsync *.[ch] Makefile $(bins) /Bossman/cz1/misc/dnapair
