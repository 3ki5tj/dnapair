remote = /Bossman/cz1/misc/dnapair

Bossman:
	rsync -avz prog/*.[Ch] prog/Makefile $(remote)/prog/
	rsync -avz doc/*.tex doc/*.pdf $(remote)/doc/
	rsync -avz NAMD/*.[Ch] NAMD/*.patch NAMD/Makefile NAMD/README.md $(remote)/NAMD/
