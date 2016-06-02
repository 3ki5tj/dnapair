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

zip: dnapairprog.zip

dnapairprog.zip::
	zip $@ *.[ch] Makefile mfrt.log mfrt.out ../doc/*.tex ../doc/*.pdf

Bossman: all
	rsync -avz *.[ch] Makefile mfrt.log $(bins) /Bossman/cz1/misc/dnapair
