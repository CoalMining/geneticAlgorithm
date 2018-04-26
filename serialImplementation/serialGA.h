#ifndef serialGA_H
#define serialGA_H

#define DEBUG_STATEMENTS 1

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class serialGA
{
private:
	string intialTempFile,finalTempFile;
	double probCrossover,probMutation,deltaX,deltaY,deltaZ,deltaT,endTime;
	int popSize,numGens;
public:	
	serialGA(){}
	serialGA(string f1,string f2):intialTempFile(f1),finalTempFile(f2){}
	void operate();
	void printOutput();
};

void serialGA::operate()
{
#if DEBUG_STATEMENTS
	cout<<"Operate function called"<<endl;
#endif
}

void serialGA::printOutput()
{
#if DEBUG_STATEMENTS
	cout<<"Print Output function called"<<endl;
#endif
}
#endif