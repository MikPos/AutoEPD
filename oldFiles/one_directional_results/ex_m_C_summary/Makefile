SHELL=/bin/bash
summary/all.tex: gen
	echo '\summaryDGHyper{\texttt{dg\_1}}{out/015_dg_1_11101}{\texttt{out/015\_dg\_1\_11101}}' >> summary/all.tex
summary.pdf: summary/all.tex
	mod_post --mode compileSummary
include summary/MakefileImpl
include summary/MakefileDep
include summary/MakefileAllTarget
include summary/MakefileClean
