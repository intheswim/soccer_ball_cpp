# requires cairo for antialiasing.
# For Ubuntu:
# sudo apt-get install libcairo2-dev

CC=g++

CFLAGS = -std=c++11 -Wall -O2 -I/usr/include/cairo/
CLIBS = -lX11 -lm -lcairo

all : main.cpp SoccerBall.h GenGraphics
		$(CC) $(CFLAGS) -o soccer main.cpp $(CLIBS) GenGraphics.o

GenGraphics: GenGraphics.cpp GenGraphics.h
		$(CC) $(CFLAGS) -c GenGraphics.cpp -o GenGraphics.o

.PHONY: clean

clean :
		-rm soccer GenGraphics.o