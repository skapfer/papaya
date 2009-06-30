unset key
plot \
	"counterexample.out/labels.out" index 0 w lp lt 1 pt 1\
	,\
	"counterexample.out/labels.out" index 1 w lp lt 1 pt 1\
	,\
	"counterexample.out/labels.out" index 2 w lp lt 1 pt 1\
	,\
	(1./0)
