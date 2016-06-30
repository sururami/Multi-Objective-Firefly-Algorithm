# Makefile for C version of I-PAES code
#
#
# To port to your machine:
#	set CFLAGS for your C compiler
#	set FLAGS for your C compiler
#	set CC to your C compiler

CC   = g++
FLAG  = -w -c -g -O
LFLAG = -lm

OBJS  = Mo_Fa_main.o paes_routines.o ran.o auxfunc.o psp_routines.o \
	print_routines.o mutation.o energy.o randomlib.o

Mo_Fa_main: Mo_Fa_main.o paes_routines.o ran.o auxfunc.o psp_routines.o \
	print_routines.o mutation.o energy.o randomlib.o

Mo_Fa_main: $(OBJS)
	$(CC) $(OBJS) $(LFLAG) -o ../Mo_Fa_main
	
Mo_Fa_main.o: global.h Mo_Fa_main.c 

	$(CC) $(FLAG) Mo_Fa_main.c
	
paes_routines.o: global.h paes_routines.c

	$(CC) $(FLAG) paes_routines.c

ran.o: ran.c

	$(CC) $(FLAG) ran.c
	
auxfunc.o: auxfunc.c

	$(CC) $(FLAG) auxfunc.c
	
psp_routines.o: psp_routines.c

	$(CC) $(FLAG) psp_routines.c
	
print_routines.o: print_routines.c

	$(CC) $(FLAG) print_routines.c
	
energy.o: energy.c

	$(CC) $(FLAG) energy.c
	
mutation.o: mutation.c

	$(CC) $(FLAG) mutation.c
	
randomlib.o: random/randomlib.c

	$(CC) $(FLAG) random/randomlib.c
	
clean:
	rm *.o ../Mo_Fa_main
