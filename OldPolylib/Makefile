CC=gcc-4.6
file=eGLJpoly.c

target:
	${CC} -c polylib.c
	${CC} -shared -fPIC -o libpolylib.dyld polylib.o 
	rm *.o
