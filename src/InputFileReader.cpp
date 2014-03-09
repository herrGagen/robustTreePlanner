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
	in_stream >> angularWidth;
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
	in_stream >> outputFilename;
	in_stream.close();

	if(!areInputsLegal() )
	{
		exit(-1);
	}
}

/**
	\brief Checks that values from input file are legal.
*/
bool InputFileReader::areInputsLegal()
{
	bool inputsAreValid = true;
	if( demandFile.size() == 0 )
	{
		std::cerr << "NULL demand file" << std::endl;	
		inputsAreValid = false;
	}
	if(weatherFileNames.size() == 0 )
	{
		std::cerr << "No weather files" << std::endl;	
		inputsAreValid = false;
	}
	for( std::vector<std::string>::iterator sIter = weatherFileNames.begin();
		sIter != weatherFileNames.end();
		++sIter)
	{
		if( sIter->size() == 0)
		{
			std::cerr << "Null weather file" << std::endl;	
			inputsAreValid = false;
		}
	}
	if( cellWidth < .0001 )
	{
		std::cerr << "Invalid cellWidth, it is < .0001" << std::endl;
		inputsAreValid = false;			
	}
	if( deviationThreshold < 0 || deviationThreshold > 1)
	{
		std::cerr << "Deviation threshold not in [0,1], it is: ";
		std::cerr << deviationThreshold << std::endl;
		inputsAreValid = false;			
	}
	if( nodeEdgeThreshold < 0 || nodeEdgeThreshold > 1)
	{
		std::cerr << "Node edge threshold not in [0,1], it is: ";
		std::cerr << nodeEdgeThreshold << std::endl;
		inputsAreValid = false;			
	}
	if( laneWidth < 0 )
	{
		std::cerr << "Warning: Lane Width is negative, so it is being read from demand files " << std::endl;
	}
	return inputsAreValid;
}