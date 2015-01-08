CC=g++
CFLAGS= -pedantic -Wall -W -O3

#CFLAGS= -pedantic -Wall -g -pg  -Wall -O0 
#// *-O -Wall


PROG1=art_1
PROG2=art_2A
PROG3=art_2A-C
PROG4=art_distance
COMMON=io_funct.o get_options.o art_common.o
OBJS1=art_1_alg.o art_1.o
OBJS2=art_2A_alg.o art_2A.o
OBJS3=art_2A-C_alg.o art_2A-C.o 
OBJS4=art_distance_alg.o art_distance.o 




art:$(COMMON) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)
	$(CC) $(CFLAGS) $(COMMON) $(OBJS1) -o $(PROG1)
	$(CC) $(CFLAGS) $(COMMON) $(OBJS2) -o $(PROG2)
	$(CC) $(CFLAGS) $(COMMON) $(OBJS3) -o $(PROG3)
	$(CC) $(CFLAGS) $(COMMON) $(OBJS4) -o $(PROG4) -lgsl -lgslcblas -lm

#static - if you choose this one - comment previous 5 lines and uncomment following 5 lines
#art:$(COMMON) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)
#	$(CC) $(CFLAGS) -static $(COMMON) $(OBJS1) -o $(PROG1) -lm -lstdc++ -lc
#	$(CC) $(CFLAGS) -static $(COMMON) $(OBJS2) -o $(PROG2) -lm -lstdc++ -lc
#	$(CC) $(CFLAGS) -static $(COMMON) $(OBJS3) -o $(PROG3) -lm -lstdc++ -lc
#	$(CC) $(CFLAGS) -static $(COMMON) $(OBJS4) -o $(PROG4) -lgsl -lgslcblas  -lm -lstdc++  -lc 



	strip $(PROG1)
	strip $(PROG2)
	strip $(PROG3)
	strip $(PROG4)

io_funct.o:io_funct.cpp io_funct.h
	$(CC) $(CFLAGS) -c io_funct.cpp -o io_funct.o

get_options.o:get_options.cpp get_options.h
	$(CC) $(CFLAGS) -c get_options.cpp -o get_options.o

art_common.o:art_common.cpp art_common.h
	$(CC) $(CFLAGS) -c art_common.cpp -o art_common.o

art_1_alg.o:art_1_alg.cpp art_1_alg.h
	$(CC) $(CFLAGS) -c art_1_alg.cpp -o art_1_alg.o

art_1.o:art_1.cpp
	$(CC) $(CFLAGS) -c art_1.cpp -o art_1.o

art_2A_alg.o:art_2A_alg.cpp art_2A_alg.h
	$(CC) $(CFLAGS) -c art_2A_alg.cpp -o art_2A_alg.o

art_2A.o:art_2A.cpp
	$(CC) $(CFLAGS) -c art_2A.cpp -o art_2A.o

art_distance_alg.o:art_distance_alg.cpp art_distance_alg.h
	$(CC) $(CFLAGS) -c art_distance_alg.cpp -o art_distance_alg.o

art_distance.o:art_distance.cpp
	$(CC) $(CFLAGS) -c art_distance.cpp -o art_distance.o

art_2A-C_alg.o:art_2A-C_alg.cpp art_2A-C_alg.h
	$(CC) $(CFLAGS) -c art_2A-C_alg.cpp -o art_2A-C_alg.o

art_2A-C.o:art_2A-C.cpp
	$(CC) $(CFLAGS) -c art_2A-C.cpp -o art_2A-C.o

clean:
	rm -f $(COMMON) $(OBJS1) $(PROG1) $(OBJS2) $(PROG2) $(OBJS3) $(PROG3) $(OBJS4) $(PROG4)
