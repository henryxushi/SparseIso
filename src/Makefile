CC=g++
CFLAGS=-O3 -I /home/cbil/Downloads/eigen-eigen-bdd17ee3b1b3 -I /home/cbil/Henry/bamtools/include -L /home/cbil/Henry/bamtools/lib -lbamtools -lboost_system -lboost_filesystem -lpthread -lboost_thread
all: SparseIso
SparseIso: SparseIso.o filter_bam.o readinstance.o utility.o Info.o options.o
		$(CC) $(CFLAGS) SparseIso.o filter_bam.o readinstance.o utility.o Info.o options.o -o SparseIso

SparseIso.o: SparseIso.cpp
						$(CC) $(CFLAGS) -c SparseIso.cpp

filter_bam.o: filter_bam.cpp filter_bam.h
						$(CC) $(CFLAGS) -c filter_bam.cpp
						
readinstance.o: readinstance.cpp readinstance.h
						$(CC) $(CFLAGS) -c readinstance.cpp

utility.o: utility.cpp utility.h
						$(CC) $(CFLAGS) -c utility.cpp				

options.o: options.cpp options.h
						$(CC) $(CFLAGS) -c options.cpp	

Info.o: Info.cpp Info.h
						$(CC) $(CFLAGS) -c Info.cpp
											
clean:
	rm -rf *.o SparseIso
