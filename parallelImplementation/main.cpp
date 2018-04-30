#include "parallelGA.h"

int main(int argc,char* argv[])
{
	string initalTempFile="";
	string finalTempFile="";

	int commSize,myRank;
	//initialize MPI env.
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

#if DEBUG_STATEMENTS
	if(myRank==0) cout<<"MPI initialized.. There are "<<commSize<<" processes."<<endl;
#endif

	serialGA sga(initalTempFile,finalTempFile,commSize,myRank);
	sga.operate();
	
	if(myRank==0) sga.printOutput();

#if DEBUG_STATEMENTS
	if(myRank==0) cout<<"Ending MPI calls.."<<endl;
#endif

	MPI_Finalize();
	return 0;
}
