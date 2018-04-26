#include "serialGA.h"

int main(int argc,char* argv[])
{
	string initalTempFile="";
	string finalTempFile="";

	if(argc<4)
	{
		cout<<"Enter the initial temperature file"<<endl;
		cin>>initalTempFile;
		cout<<"Enter the final temperature file"<<endl;
		cin>>finalTempFile;
	}else
	{
		initalTempFile = argv[1];
		finalTempFile = argv[0];
	}

	serialGA sga(initalTempFile,finalTempFile);
	sga.operate();
	sga.printOutput();

	return 0;
}
