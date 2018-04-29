#include "serialGA.h"

int main(int argc,char* argv[])
{
	string initalTempFile="";
	string finalTempFile="";

	serialGA sga(initalTempFile,finalTempFile);
	sga.operate();
	sga.printOutput();

	return 0;
}
