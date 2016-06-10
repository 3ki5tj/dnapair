src := $(wildcard *.h *.C)

srcdir = ~/lwork/NAMD_2.11_Source/src
destdir = ~/lwork/NAMD_mods/NAMD_2.11_dnapair/src

NAMD.patch::
	$(MAKE) ln
	diff -u -wdB --exclude=".*" --exclude="*~" $(srcdir) $(destdir) > $@

# make links
ln::
	ln -f `diff -qrwdB --exclude=".*" --exclude="*~" $(srcdir) $(destdir) | cut -f 4  -d" "` .

Bossman: NAMD.patch
	$(MAKE) -C .. $@

clean:
	rm -rf *~ .*~
