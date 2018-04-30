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
	popSize = 4000;
	numGens = 5000;
	endTime = 1.0;

	dim1 = 11;
	dim2 = 11;

	initialTempData.resize(dim2*dim1,0);
	refTempData.resize(dim2*dim1,0);
	tempTempData.resize(dim2*dim1,0);
	finalTempData.resize(dim2*dim1,0);

	bestFinalTempApprox.resize(dim2*dim1,0);

	population.resize(popSize*dim1*dim2,0);
	nextPopulation.resize(popSize*dim1*dim2,0);
	objFunction.resize(popSize,0);
	nextObjFunction.resize(popSize,0);
	bestMember.resize(dim1*dim2,0);
	bestMemberRMSE.resize(numGens,0);

	bestObjFunction = std::numeric_limits<double>::max();
	bestObjFunctions.resize(numProcesses);
	generate(bestObjFunctions.begin(),bestObjFunctions.end(),std::numeric_limits<double>::max());

#if DEBUG_STATEMENTS
	cout<<"Defining the parameters successful.."<<endl;
#endif
}

void serialGA::saveMatrix(const vector<double> &source, vector<double> &destination)
{
	for(int i=0;i<dim1*dim2;i++)
	{
		destination[i] = source[i];
	}
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

	saveMatrix(finalTempData,bestFinalTempApprox);

#if DEBUG_STATEMENTS
	cout<<"Printing the best matix being passed from one generation to other"<<endl;
	printMatrix(bestMember);
#endif
}

//this function calculates the final temperature over the time
// time ranges from t to t+endTime with an increase of deltaT step
// the result will be filled in the finalTempData
// k is the population that we are using now
void serialGA::calculateTemp(double t,int index)
{
	double finalTime = t+endTime;
	int startIndex = index*dim1*dim2;

	if(t==0.0)
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
	sumDiff = sqrt(sumDiff)/(dim2*dim1);
#if DEBUG_STATEMENTS
	cout<<"fitness Value is: "<<sumDiff<<" ("<<bestObjFunction<<")"<<endl;
#endif
	return sumDiff;
}

void serialGA::saveInitPopMember(int myIndex)
{
	int startIndex = myIndex*dim1*dim2;

	for(int i=0;i<dim1*dim2;i++)
	{
		bestMember[i] = population[startIndex+i];
	}
	bestObjFunction = objFunction[myIndex];
}

void serialGA::initPopulation()
{
	//use generate ti fill the elements in the population with random numbers
	std::srand(time(0));
	std::generate(population.begin(),population.end(),generator);

#if DEBUG_STATEMENTS
	cout<<"Printing one random instance of population"<<endl;
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
		if(objFunction[i]<bestObjFunction)
		{
			saveInitPopMember(i);
		}
	}
#if DEBUG_STATEMENTS
	cout<<"The initial best matrix is:"<<endl;
	printMatrix(bestMember);
#endif
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

	//read at process 0 and broadcast to all other processes
	if(myProcessRank==0)
	{
		//make sure that the file names are there
		if(initialTempFile==""||refTempFile=="")
		{
			cout<<"The file names are empty...Please provide "<<endl;
			initialTempFile = "U_InitialTemp.dat";
			cout<<"\n:Initial Temp file name is "<<initialTempFile<<endl;
			refTempFile = "U_FinalTemp.dat";
			cout<<"\n:Reference Temp file name is "<<refTempFile<<endl;
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
		cout<<"The read refTempData matrix is::"<<endl;
		printMatrix(refTempData);
#endif

		//initialize the population
		// with random values for the first generation, other generations will derive from this generaiton
		initPopulation();
	}

	//send the initial and final temp data to all the processes
	MPI_Bcast(&initialTempData,initialTempData.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&finalTempData,initialTempData.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);

	//send the initial population to other processes
	MPI_Bcast(&population,population.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	//send the objective function to other processes
	MPI_Bcast(&objFunction,objFunction.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&bestObjFunction,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#if DEBUG_STATEMENTS
	char checkInput = 0;
	srand(time(0));
	int procTemp1 = rand()%numProcesses;
	int procTemp2 = rand()%numProcesses;
	int popIndexTemp = rand()%popSize;
	if(myProcessRank==0)
	{
		cout<<"Broadcast to other processes completed"<<endl;
		cout<<"Checking one random instance "<< popIndexTemp<<" of population two random processes"<< procTemp1<<" and "<<procTemp2<<" to verify"<<endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myProcessRank==procTemp1)
	{
		cout<<"Printing population from rank..:"<<myProcessRank<<endl;
		for(int i=0;i<dim1;i++)
		{
			for(int j=0;j<dim2;j++)
			{
				cout<<population[(popIndexTemp*dim1*dim2)+(i*dim2+j)]<<" ";
			}
			cout<<endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myProcessRank==procTemp2)
	{
		cout<<"Printing population from rank..:"<<myProcessRank<<endl;
		for(int i=0;i<dim1;i++)
		{
			for(int j=0;j<dim2;j++)
			{
				cout<<population[(popIndexTemp*dim1*dim2)+(i*dim2+j)]<<" ";
			}
			cout<<endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	checkInput = getchar();
#endif

	int parentIndex1,parentIndex2;
	int populationPerProcess = (popSize/numProcesses);
	int myPopulationStartIndex = myProcessRank*populationPerProcess;
	int myPopulationEndIndex = (myProcessRank+1)*populationPerProcess;		

	//perform main genetic algorithm operations here
	//the operation will go for numGens times so..
	for(int i=0;i<numGens;i++)
	{

#if DEBUG_STATEMENTS
		if(myProcessRank==0)
		{
			cout<<"*********************************************"<<endl;
			cout<<"*********Generation "<<i<<" started**********"<<endl;
			cout<<"*********************************************"<<endl;
		}
#endif
		//perform for each member of the current generation ie for each population
		//for all the other members of the population, perform GA
		for(int j=myPopulationStartIndex;j<myPopulationEndIndex;j++)
		{
			if(j==0)
			{
				//for the first member of this generation, let us copy the best member from the previous generation
				//this method is called elitism
				copyPrevBest(0);
			}else{
				//first select best parents
				selectParents(parentIndex1,parentIndex2,j);
				crossOver(parentIndex1,parentIndex2,j);
				mutate(j);
				calculateTemp(0.0,j);

				nextObjFunction[j] = fitnessFunction();
				if(nextObjFunction[j]<bestObjFunction)
				{
#if DEBUG_STATEMENTS
					cout<<"The best member is changed...index is "<<j<<" and fitness value is "<<nextObjFunction[j]<<endl;
#endif
					savePopMember(j);
				}
			}
		}

		bestObjFunctions[myProcessRank] = bestObjFunction;

		//broadcast the population data to all other processes
		MPI_Barrier(MPI_COMM_WORLD);
		for(int idx = 0;idx<numProcesses;idx++)
		{
			int bcastStartIndex = idx*dim1*dim2*populationPerProcess;
			MPI_Bcast(&nextPopulation[bcastStartIndex],populationPerProcess*dim1*dim2,MPI_DOUBLE,idx,MPI_COMM_WORLD);	
			MPI_Bcast(&nextObjFunction[populationPerProcess*idx],populationPerProcess,MPI_DOUBLE,idx,MPI_COMM_WORLD);	
			MPI_Bcast(&bestObjFunctions[idx],1,MPI_DOUBLE,idx,MPI_COMM_WORLD);	
#if DEBUG_STATEMENTS
			if(myProcessRank==idx)
			{
				cout<<"Broadcasting next population and next object function from process"<< myProcessRank<<" to all other processes completed"<<endl;
			}
#endif	
		}

		//in rank=0, find the least of the objective functions and send again to all other processes
		int minPosition;
		if(myProcessRank==0)
		{
			minPosition = min_element(bestObjFunctions.begin(),bestObjFunctions.end())-bestObjFunctions.begin();
		}
		MPI_Bcast(&minPosition,1,MPI_INT,0,MPI_COMM_WORLD);	
		bestObjFunction = bestObjFunctions[minPosition];		
		MPI_Bcast(&bestMember,dim1*dim2,MPI_DOUBLE,minPosition,MPI_COMM_WORLD);
		//copy the population and objective function to next population
		nextObjFunction.swap(objFunction);
		nextPopulation.swap(population);

		MPI_Barrier(MPI_COMM_WORLD);

	}
#if DEBUG_STATEMENTS
	if(myProcessRank==0)
	{
		cout<<"The final best approximation of temperature is "<<endl;
		printMatrix(bestFinalTempApprox);
	}
#endif
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
#if DEBUG_STATEMENTS
	cout<<"Writing best RMSE for each generation to a file"<<endl;
	fstream bestRMSE("bestRMSE.csv",ios::out);
	for(int i=0;i<numGens;i++)
	{
		bestRMSE<<i<<","<<bestMemberRMSE[i]<<endl;
	}
	bestRMSE.close();
	cout<<"Writing completed"<<endl;
#endif
}

void serialGA::savePopMember(int myIndex)
{
	int startIndex = myIndex*dim1*dim2;

	for(int i=0;i<dim1*dim2;i++)
	{
		bestMember[i] = nextPopulation[startIndex+i];
	}
	bestObjFunctions[myProcessRank] = nextObjFunction[myIndex];
}

void serialGA::mutate(int myIndex)
{
	//rand in range 0 to 1 double
	int myNum = rand();
	double nextNum = (double) myNum/std::numeric_limits<int>::max();

	int mutationPoint;

	if(nextNum<probMutation)
	{
#if DEBUG_STATEMENTS
		cout<<"Mutation occurred"<<endl;
#endif
		mutationPoint = rand()%popSize;
		nextPopulation[myIndex*dim1*dim2] = generator();
	}
}

void serialGA::crossOver(int parentIndex1,int parentIndex2,int myIndex)
{
	int startIndex = myIndex*dim1*dim2;
	int crossPoint;

//	random number in range 0 to 1 double
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

#if DEBUG_STATEMENTS
	cout<<"Parents selected are: "<<parentIndex1<<" and "<<parentIndex2<<endl;
#endif
}