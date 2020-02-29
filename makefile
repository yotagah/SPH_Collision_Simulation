all:
	gcc *.c -o sph -lm -Wall

clean:
	-rm -f *~ *.o

purge: clean
	-rm -f sph