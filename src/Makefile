CC=g++
PATH_EIGEN=/home/cbil/Downloads/eigen-eigen-bdd17ee3b1b3
PATH_BAMTOOLS_INCLUDE=/home/cbil/Henry/bamtools/include
PATH_BAMTOOLS_LIB=/home/cbil/Henry/bamtools/lib

CFLAGS=-O3 -I $(PATH_EIGEN) -I $(PATH_BAMTOOLS_INCLUDE) -L $(PATH_BAMTOOLS_LIB) -lbamtools -lboost_system -lboost_filesystem -lpthread -lboost_thread
all: SparseIso
SparseIso: SparseIso.o filter_bam.o readinstance.o utility.o Info.o options.o getinsertsize.o
		$(CC) $(CFLAGS) SparseIso.o filter_bam.o readinstance.o utility.o Info.o options.o getinsertsize.o -o ../SparseIso

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
						
getinsertsize.o: getinsertsize.cpp getinsertsize.h
						$(CC) $(CFLAGS) -c getinsertsize.cpp						
											
clean:
	rm -rf *.o SparseIso
