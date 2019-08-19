#important: cuba must have been compiled with -fPIC option

all:
	gcc -O3 -fPIC -shared -o c_speedup.so c_speedup.c -llapack -lblas
	#gcc -O3 -fPIC -shared -o c_speedup.so c_speedup.c -lmkl_core -lmkl_intel_lp64

