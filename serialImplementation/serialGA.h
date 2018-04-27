#ifndef serialGA_H
#define serialGA_H

#define DEBUG_STATEMENTS 1

#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;


//this is the generator function to generate random double number
double generator()
{
	int myNum = std::rand();
	return (double)myNum* 0.1/std::numeric_limits<int>::max();
}

class serialGA
{
private:
	string initialTempFile,refTempFile;
	
	double probCrossover,probMutation,deltaX,deltaY,deltaT,endTime,bestObjFunction;
	int popSize,numGens,dim1,dim2;

	vector<double> initialTempData,refTempData,tempTempData,finalTempData;
	vector<double> population,objFunction,nextPopulation,nextObjFunction;
	vector<double> bestMember;
public:	
	serialGA()
	{
		defineParameters();
	}
	serialGA(string f1,string f2):initialTempFile(f1),refTempFile(f2)
	{
		defineParameters();
	}
	~serialGA(){}
	void operate();
	bool readMatrix(vector<double> &matrix,string fromFile);
	void printOutput();
	void defineParameters();
	void copyPrevBest(int myIndex);
	void calculateTemp(double t,int index);
	double fitnessFunction();
	void initPopulation();
	void savePopMember(int myIndex);
	void mutate(int myIndex);
	void crossOver(int parentIndex1,int parentIndex2,int myIndex);
	void selectParent(int &indexParent);
	void selectParents(int &parentIndex1, int &parentIndex2, int indexMine);
};

#endif