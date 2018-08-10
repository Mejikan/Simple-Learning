all: learn.o matrix.o
	gcc -Wall -Werror -fsanitize=address -o learn learn.o matrix.o

learn.o: learn.c matrix.h
	gcc -c learn.c

matrix.o: matrix.c matrix.h
	gcc -c matrix.c

clean:
	rm -f learn learn.o matrix.o
