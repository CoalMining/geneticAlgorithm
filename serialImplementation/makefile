#########################################################3
###############3 This is the makefile for the serial project
############ g++ main.cpp serialGA.h -o serial ###########
##########################################################
CC=g++
DEPS=serialGA.h
CFLAGS=
OBJ=serialGA.o main.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

serial: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)