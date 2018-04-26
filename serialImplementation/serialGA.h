#ifndef serialGA_H
#define serialGA_H

#define DEBUG_STATEMENTS 1

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <ctime>

using namespace std;

class serialGA
{
private:
	string initialTempFile,finalTempFile;
	
	double probCrossover,probMutation,deltaX,deltaY,deltaZ,endTime,bestObjFunction;
	int popSize,numGens,dim1,dim2;

	vector<double> initialTempData,finalTempData;
	vector<double> population,objFunction,nextPopulation,nextObjFunction;
	vector<double> bestMember;
public:	
	serialGA()
	{
		defineParameters();
	}
	serialGA(string f1,string f2):initialTempFile(f1),finalTempFile(f2)
	{
		defineParameters();
	}
	~serialGA(){}
	void operate();
	bool readMatrix(vector<double> &matrix,int d1,int d2,string fromFile);
	void printOutput();
	void defineParameters();
};

//this function defines the initial parameters that are to be constant
void serialGA::defineParameters()
{
#if DEBUG_STATEMENTS
	cout<<"Defining the parameters that will remain constant throughout...."<<endl;
#endif
	probCrossover = 0.90;
	probMutation = 0.05;
	deltaX = 0.1;
	deltaY = 0.1;
	deltaT = 0.01;
	popSize = 400;
	numGens = 2000;
	endTime = 1.0;

	dim1 = 11;
	dim2 = 11;

	initialTempData.reserve(dim2*dim1);
	finalTempData.reserve(dim2*dim1);

	population.reserve(popSize*dim1*dim2);
	nextPopulation.reserve(popSize*dim1*dim2);
	objFunction.reserve(popSize);
	nextObjFunction.reserve(popSize);
	bestMember.reserve(dim1*dim2);

#if DEBUG_STATEMENTS
	cout<<"Defining the parameters successful.."<<endl;
#endif
}

//read fron file fromFile to the matrix
// number of items read will be d1*d2
bool serialGA::readMatrix(vector<double> &matrix,int& d1,int& d2,string fromFile)
{
	fstream myFile(fromFile.c_str());
	if(!myFile.is_open())
	{
		cout<<"Error in reading the file "<<fromFile<<endl;
		return false;
	}
#if DEBUG_STATEMENTS
	cout<<"Reading from the file "<<fromFile<<endl;
#endif

	double tempVal  = 0.0;
	for(int i=0;i<d1;i++)
	{
		try
		{
			for(int j=0;j<d2;j++)
			{
				myFile>>tempVal; 
				matrix.push_back(tempVal);
			}
		}
		catch(...)
		{
			cout<<"The error might be because the matrix is not properly sized..\n\t while reading "<<fromFile<<endl;
			return false;
		}
	}
#if DEBUG_STATEMENTS
	cout<<"Reading from the file "<<fromFile<<"completed"<<endl;
#endif
	return true;
}

//this function copies the population that has best objFunction value from the prevGen to firstItem of nextGen
//copyPosition is the position or the population index that will be the copy of best member
void serialGA::copyPrevBest(vector<double> &child,int copyPosition,int d1,int d2)
{
	for(int i=0;i<d1*d2;i++)
	{
		child[copyPosition*d1*d2+i] = bestMember[i];
	}
	nextObjFunction[copyPosition] = bestObjFunction;
}
//this is the main function that operates the GA 
void serialGA::operate()
{
#if DEBUG_STATEMENTS
	cout<<"Operate function called"<<endl;
#endif

	//make sure that the file names are there
	if(initialTempFile==""||finalTempFile=="")
	{
		cout<<"The file names are empty...Please provide "<<endl;
		cout<<"\n:Initial Temp file name"<<endl;
		cin>>initialTempFile;
		cout<<"\n:Final Temp file name"<<endl;
		cin>>finalTempFile;
	}

	//read the data from files
	if(!readMatrix(initialTempData,dim1,dim2,initialTempFile))
	{
		cout<<"Initial temperature matrix read unsuccessful...."<<endl;
		return;
	}
	if(!readMatrix(finalTempData,dim1,dim2,finalTempFile))
	{
		cout<<"Final temperature matrix read unsuccessful...."<<endl;
		return;
	}

	//perform main genetic algorithm operations here
	//the operation will go for numGens times so..
	for(int i=0;i<numGens;i++)
	{
		//perform for each member of the current generation ie for each population
		for(int j=0;j<popSize;j++)
		{
			//for the first member of this generation, let us copy the best member from the previous generation
			if(j==0)
			{
				//this method is called elitism
				copyPrevBest(nextPopulation,j,dim1,dim2);
			}else
			{

			}
		}
	}
}

//print the final matrix output
void serialGA::printOutput()
{
#if DEBUG_STATEMENTS
	cout<<"Print Output function called"<<endl;
#endif
}
#endif