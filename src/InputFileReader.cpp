#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include "InputFileReader.h"

/**
   \brief Simple constructor that parses an input file
   
   To see the format of the input file, see the class description.
   
   \param inputFileName Name of input file (i.e. input.txt)
*/
InputFileReader::InputFileReader(std::string inputFileName)
{
	std::ifstream in_stream;
	std::string line;

	std::cout << "Input file: " << inputFileName.c_str() << std::endl;
	in_stream.open(inputFileName.c_str());

	in_stream >> demandFile;
	
	unsigned int numDataFiles;
	in_stream >> numDataFiles;
	weatherFileNames.resize(numDataFiles);
	for(unsigned int i = 0; i<numDataFiles; i++)
	{
		std::string temp;
		in_stream >> temp;
		weatherFileNames[i] = temp;
	}

	in_stream >> cellWidth;
	in_stream >> quadrantAngle;
	in_stream >> deviationThreshold;
	in_stream >> nodeEdgeThreshold;
	in_stream >> angleOffset;
	in_stream >> laneWidth;
	in_stream >> numFixedNodes;
	in_stream >> operFlex1;
	in_stream >> operFlex2;
	in_stream >> operFlex3;
	in_stream >> timestamp;
	in_stream.close();
}