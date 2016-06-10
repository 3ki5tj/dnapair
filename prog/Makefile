GP = $(wildcard *.gp)
EPS = $(patsubst %.gp,%.eps, $(GP))
PDF = $(patsubst %.gp,%.pdf, $(GP))

CC = icc -Wall -Wremarks -O2 -g
LM =
#CC = gcc -O3 -Wall -Wextra -g
#LM = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
deps = $(wildcard *.h) Makefile

all: $(bins) $(PDF)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

clean:
	rm -rf $(bins) \
	  a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos
	rstrip.py -Rlv

$(PDF) : %.pdf : %.eps $(epsdeps)
	epstopdf $<

$(EPS) : %.eps : %.gp $(epsdeps)
	gnuplot $<

zip: dnapairprog.zip

dnapairprog.zip::
	$(MAKE) $(PDF)
	zip $@ *.[ch] Makefile mfrt.log mfrt.out $(GP) $(EPS) $(PDF) ../doc/*.tex ../doc/*.pdf

Bossman: all
	rsync -avz *.[ch] Makefile mfrt.log $(bins) $(GP) $(EPS) $(PDF) /Bossman/cz1/misc/dnapair
