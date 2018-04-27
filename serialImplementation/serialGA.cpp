#include "serialGA.h"

//this is the generator function to generate random double number
double generator()
{
	int myNum = std::rand();
	return (double)myNum* 0.1/std::numeric_limits<int>::max();
}

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

	initialTempData.resize(dim2*dim1,0);
	refTempData.resize(dim2*dim1,0);
	tempTempData.resize(dim2*dim1,0);
	finalTempData.resize(dim2*dim1,0);

	population.resize(popSize*dim1*dim2,0);
	nextPopulation.resize(popSize*dim1*dim2,0);
	objFunction.resize(popSize,0);
	nextObjFunction.resize(popSize,0);
	bestMember.resize(dim1*dim2,0);

#if DEBUG_STATEMENTS
	cout<<"Defining the parameters successful.."<<endl;
#endif
}


//read fron file fromFile to the matrix
// number of items read will be d1*d2
bool serialGA::readMatrix(vector<double> &matrix,string fromFile)
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
	for(int i=0;i<dim1;i++)
	{
		try
		{
			for(int j=0;j<dim2;j++)
			{
				myFile>>tempVal; 
				matrix[i*dim2+j] = tempVal;
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

//this function copies the best member to the next population member
// position to be copied in the next population is identified by copyPosition
void serialGA::copyPrevBest(int myIndex)
{
	int startIndex = myIndex*dim1*dim2;
	for(int i=0;i<dim1*dim2;i++)
	{
		nextPopulation[startIndex+i] = bestMember[i];
	}
	nextObjFunction[myIndex] = bestObjFunction;
}

//this function calculates the final temperature over the time
// time ranges from t to t+endTime with an increase of deltaT step
// the result will be filled in the finalTempData
// k is the population that we are using now
void serialGA::calculateTemp(double t,int index)
{
	double finalTime = t+endTime;
	int startIndex = index*dim1*dim2;

	if(t==0)
	{
		//for the initial time instance
		for(int i=0;i<dim1;i++)
		{
			for(int j=0;j<dim2;j++)
			{
				tempTempData[i*dim2+j] = initialTempData[i*dim2+j];
			}
		}
	}
	for(;t<=finalTime;t+=deltaT)
	{
		for(int i=1;i<dim1-1;i++)
		{
			for(int j=1;j<dim2-1;j++)
			{
				finalTempData[i*dim2+j] = tempTempData[i*dim2+j]+population[startIndex+i*dim2+j]*deltaT*
				(
					(
						(tempTempData[(i+1)*dim2+j]-2.0*tempTempData[i*dim2+j]+tempTempData[(i-1)*dim2+j])/(deltaX*deltaX)
					)+
					(
						(tempTempData[(i)*dim2+(j+1)]-2.0*tempTempData[i*dim2+j]+tempTempData[(i)*dim2+(j-1)])/(deltaY*deltaY)
					)
				);
			}
		}

		for(int i=1;i<dim1-1;i++)
		{
			for(int j=1;j<dim2-1;j++)
			{
				tempTempData[i*dim2+j] = finalTempData[i*dim2+j];
			}
		}
	}
}

double serialGA::fitnessFunction()
{
	double sumDiff = 0.0;
	for(int i=0;i<dim2*dim1;i++)
	{
		double diff = refTempData[i]-finalTempData[i];
		sumDiff+=(diff*diff);
	}
	return sqrt(sumDiff)/(dim2*dim1);
}

void serialGA::initPopulation()
{
	//use generate ti fill the elements in the population with random numbers
	std::srand(time(0));
	std::generate(population.begin(),population.end(),generator);

#if DEBUG_STATEMENTS
	cout<<"Printing one randon instance of population"<<endl;
	int idx = rand()%popSize;
	for(int i=0;i<dim1;i++)
	{
		for(int j=0;j<dim2;j++)
		{
			cout<<population[idx*dim1*dim2+(i*dim2+j)]<<" ";
		}
		cout<<endl;
	}
#endif

	for(int i=0;i<popSize;i++)
	{
		double timeStart = 0.0;
		//now calcualte the population or temperature
		//this function goes through a series of iterations over time intervals and calculates
		//.. the final temperature that will be for next generation
		calculateTemp(0.0,i);
		//this fitness functin is based on the tempTempData and refTempData
		objFunction[i] = fitnessFunction();
	}
}

void serialGA::printMatrix(vector<double> &myMatrix)
{
	for(int i=0;i<dim1;i++)
	{
		for(int j=0;j<dim2;j++)
		{
			cout<<myMatrix[i*dim2+j]<<" ";
		}
		cout<<endl;
	}
}

//this is the main function that operates the GA 
void serialGA::operate()
{
#if DEBUG_STATEMENTS
	cout<<"Operate function called"<<endl;
#endif

	//make sure that the file names are there
	if(initialTempFile==""||refTempFile=="")
	{
		cout<<"The file names are empty...Please provide "<<endl;
		cout<<"\n:Initial Temp file name"<<endl;
		cin>>initialTempFile;
		cout<<"\n:Reference Temp file name"<<endl;
		cin>>refTempFile;
	}

	//read the data from files
	if(!readMatrix(initialTempData,initialTempFile))
	{
		cout<<"Initial temperature matrix read unsuccessful...."<<endl;
		return;
	}
#if DEBUG_STATEMENTS
	cout<<"The read initialTempData matrix is::"<<endl;
	printMatrix(initialTempData);
#endif
	if(!readMatrix(refTempData,refTempFile))
	{
		cout<<"Reference temperature matrix read unsuccessful...."<<endl;
		return;
	}
#if DEBUG_STATEMENTS
	cout<<"The read finalTempData matrix is::"<<endl;
	printMatrix(finalTempData);
#endif

	//initialize the population
	// with random values for the first generation, other generations will derive from this generaiton
	initPopulation();

	int parentIndex1,parentIndex2;
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
				copyPrevBest(j);
			}else
			{
				//for all the other members of the population, perform GA

				//first select best parents
				selectParents(parentIndex1,parentIndex2,j);
				crossOver(parentIndex1,parentIndex2,j);
				mutate(j);
				calculateTemp(0.0,j);

				nextObjFunction[j] = fitnessFunction();
				if(nextObjFunction[j]>bestObjFunction)
				{
					savePopMember(j);
				}
			}
		}

		//copy the population and objective function to next population
		nextObjFunction.swap(objFunction);
		nextPopulation.swap(population);
	}
}

void serialGA::printOutput()
{
#if DEBUG_STATEMENTS
	cout<<"Print Output function called"<<endl;
#endif
	cout<<"The final obtained K matrix is "<<endl;
	for(int i=0;i<dim1;i++)
	{
		for(int j=0;j<dim2;j++)
		{
			cout<<bestMember[i*dim2+j]<<" ";
		}
		cout<<endl;
	}
}

void serialGA::savePopMember(int myIndex)
{
	int startIndex = myIndex*dim1*dim2;

	for(int i=0;i<dim1*dim2;i++)
	{
		bestMember[i] = nextPopulation[startIndex+i];
	}
	bestObjFunction = nextObjFunction[myIndex];
}

void serialGA::mutate(int myIndex)
{
	int myNum = rand();
	double nextNum = (double) myNum/std::numeric_limits<int>::max();

	int mutationPoint;

	if(nextNum<probMutation)
	{
		mutationPoint = rand()%popSize;
		nextPopulation[myIndex*dim1*dim2] = generator();
	}
}

void serialGA::crossOver(int parentIndex1,int parentIndex2,int myIndex)
{
	int startIndex = myIndex*dim1*dim2;
	int crossPoint;

	int myNum = rand();
	double nextNum = (double) myNum/std::numeric_limits<int>::max();

	if(nextNum<probCrossover)
	{
		//then crossover occurs
		crossPoint = rand()%popSize;
		for(int i=0;i<dim1*dim2;i++)
		{
			if(i<crossPoint)
			{
				nextPopulation[startIndex+i] = population[parentIndex1*dim1*dim2+i];
			}else
			{
				nextPopulation[startIndex+i] = population[parentIndex2*dim1*dim2+i];
			}
		}

	}else
	{
		//no crossover, just copy from one parent
		for(int i=0;i<dim1*dim2;i++)
		{
			nextPopulation[startIndex+i] = population[parentIndex1*dim1*dim2+i];
		}
	}
}

void serialGA::selectParent(int &indexParent)
{
	int p1,p2;
	
	p1 = rand()%popSize;
	p2 = rand()%popSize;
	if(objFunction[p1]>objFunction[p2])
		indexParent = p1;
	else
		indexParent = p2;
}
void serialGA::selectParents(int &parentIndex1, int &parentIndex2, int indexMine)
{
	selectParent(parentIndex1);
	do
	{
		selectParent(parentIndex2);
	}while(parentIndex1==parentIndex2);
}