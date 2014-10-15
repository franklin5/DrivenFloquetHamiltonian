CC=mpicxx
eigenbin = /home/ld7/bin	
CFLAGS=-c -Wall -I${eigenbin}
LFLAGS=-limf -lm
all: floquet
#entry.o floquet.o lgwt.o ioFloquet.o distribution.o BdGmatrix.o
floquet: entry.o floquet.o lgwt.o ioFloquet.o distribution.o BdGmatrix.o
	$(CC) *.o -o floquet $(LFLAGS)

entry.o: entry.cpp
	$(CC) $(CFLAGS) entry.cpp

floquet.o: floquet.cpp
	$(CC) $(CFLAGS) floquet.cpp

ioFloquet.o: ioFloquet.cpp
	$(CC) $(CFLAGS) ioFloquet.cpp

distribution.o: distribution.cpp
	$(CC) $(CFLAGS) distribution.cpp

BdGmatrix.o: BdGmatrix.cpp
	$(CC) $(CFLAGS) BdGmatrix.cpp

lgwt.o: lgwt.cpp
	$(CC) $(CFLAGS) lgwt.cpp

touch: 
	touch *.cpp *.h
clean:
	rm *.o floquet *~ *#

