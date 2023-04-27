CC = gcc
LIB = -lm -lgsl -lgslcblas
DEBUG = -g

test : *.c
	$(CC) -o test *.c $(LIB)
debug : *.c
	$(CC) -o test *.c $(LIB) $(DEBUG)
clean :
	rm *.o
