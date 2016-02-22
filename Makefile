CC=g++

all:
	$(CC) -Wall -std=c++11 -O2 -fPIC -c cgb.cpp -I/usr/include/python2.7
	$(CC) -shared -o cgb.so cgb.o
