#########################################################3
###############3 This is the makefile for the serial project
############ g++ main.cpp serialGA.h -o serial ###########
##########################################################
CC=g++
MPICC=mpicc
DEPS=parallelGA.h
CFLAGS=
OBJ=parallelGA.o main.o

%.o: %.cpp $(DEPS)
	$(MPICC) -c -o $@ $< $(CFLAGS)

serial: $(OBJ)
	$(MPICC) -o $@ $^ $(CFLAGS)