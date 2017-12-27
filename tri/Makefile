#This Makefile compiles
#tri.c (the default target) and comb.c
#
#-----------------------------------------
tri:	tri.o
	gcc -o tri tri.o 
tri.o:  tri.c macros.h
	gcc -c tri.c

#-----------------------------------------
comb:	comb.o
	gcc -o comb comb.o 
comb.o:  comb.c
	gcc -c comb.c
#-----------------------------------------
