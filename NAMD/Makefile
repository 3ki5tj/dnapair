src := $(wildcard *.h *.C)

NAMD.patch::
	echo $(src)
	diff -u --exclude=".*" --exclude="*~" \
	  ~/lwork/NAMD_2.11_Source/src \
	  ~/lwork/NAMD_mods/NAMD_2.11_dnapair/src > $@

Bossman:
	$(MAKE) -C .. $@

clean:
	rm -rf *~ .*~
