#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

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
	if(!in_stream.is_open() )
	{
		std::cerr << "Input file " << inputFileName << " not on path.  Aborting." << std::endl;
		std::exit(-1);
	}

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
	in_stream >> demandShift;
	in_stream >> demandDrop;
	unsigned int numOperFlex;
	in_stream >> numOperFlex;
	operFlex.resize(numOperFlex);
	for(unsigned int i = 0; i < numOperFlex; i++)
	{
		in_stream >> operFlex[i];
	}
	std::sort(operFlex.begin(), operFlex.end() );
	in_stream >> timestamp;
	in_stream.close();
}
