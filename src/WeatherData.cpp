#include "WeatherData.h"
#include "math.h"
#include <algorithm>
#include "RoutingDAG.h"
#include <time.h>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>

WeatherData::WeatherData(void)
{
	numPoints = 0;
	minAlt = maxAlt = 0;				
	minProbDev = maxProbDev =0;	
	cellWidth = 0;
	cellHeight = 0;
	probability = 0;
}

WeatherData::~WeatherData(void)
{
	xCoors.clear();
	yCoors.clear();
	altitudes.clear();
	probDeviation.clear();
}

void WeatherData::reset()		//reset the weather data
{
	numPoints = 0;
	minAlt = maxAlt = 0;				
	minProbDev = maxProbDev =0;	
	xCoors.clear();
	yCoors.clear();
	altitudes.clear();
	probDeviation.clear();
	cellWidth = 0;
	cellHeight = 0;
	probability = 0;
}

// the content of the weather data file is stored in the char array buffer
// this function analyses the data and store the data members into their vectors
// the range of a box is given so that anything that is not in this range will NOT be stored in the vectors (weather trimming)
bool WeatherData::readInFileData(std::string fileName, double rangeMinLati, double rangeMinLong, double rangeMaxLati, double rangeMaxLong)
{
	reset();											// first reset the vectors to empty, then read in

	const char *fname = fileName.c_str();
	std::ifstream dataStream(fname, std::ifstream::in );
	size_t fileSize = 0;
	dataStream.seekg(0, std::ios_base::end);
	fileSize = (size_t)dataStream.tellg();
	dataStream.seekg(0, std::ios_base::beg);

	// Due to incompatability across multiple platforms, this section is being 
	// redone with more modern c++ tools than pointers, and strcomp()

	std::string thisLine;
	for(unsigned int i = 0; i<5; i++) // Try and find the "Probability" line in the first 5 lines of the file.
	{
		std::getline(dataStream, thisLine);
		size_t found = thisLine.find("Probability");
		if(found != std::string::npos)
		{
			unsigned beginNumber = thisLine.find_first_of("-.0123456789",found);
			unsigned endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);
			endNumber = (endNumber < thisLine.length() ) ? endNumber : thisLine.length()-1;
			std::string tempString = thisLine.substr(beginNumber, endNumber);
			probability = (double) ::atof( (thisLine.substr(beginNumber, endNumber-beginNumber)).c_str() );
			break; // leave the for loop once probability has been found.
		}
		if(i==4)
		{
			std::cerr << "The input file had no probability associated with it." << std::endl;
			return( handleInputData() );
		}
	}

	// New format (2013-03-14):
	// Each line of this file after this point contains 7 numbers:
	// someNum, anotherNum, anotherNum, xCoors, yCoors, altitudes, probDeviation
	// The Coordinates have to be checked with:
	// (tempX>rangeMinLati && tempX<rangeMaxLati && tempY>rangeMinLong && tempY<rangeMaxLong)
	// to ensure that we should push this element into the arrays.
	while(!dataStream.eof() )
	{
		std::getline(dataStream, thisLine);
		double values[4]; // Storage for tempX, tempY, tempAltitude, tempProbability		
		size_t beginNumber = thisLine.find_first_of("-.0123456789",0);
		// Just a double check to ensure we aren't feeding the parser garbage lines.
		if(beginNumber == std::string::npos)
		{
			break;
		}
		size_t endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);

		// Ignore the first three values of each line, since they are not used by this program.
		for(unsigned int i = 0; i < 3; i++)
		{
			beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
			endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);
		}

		for(unsigned int i = 0; i<4; i++)
		{
			std::string tempString = thisLine.substr(beginNumber, endNumber-beginNumber);
			values[i] = (double) ::atof( tempString.c_str() );
			// Move the search window for a number in the current line.
			beginNumber = thisLine.find_first_of("-.0123456789",endNumber);
			endNumber = thisLine.find_first_not_of("-.0123456789",beginNumber);
		}

		if(values[0]>rangeMinLati && 
			values[0]<rangeMaxLati && 
			values[1]>rangeMinLong && 
			values[1]<rangeMaxLong)
		{
			xCoors.push_back(values[0]);	
			yCoors.push_back(values[1]);
			altitudes.push_back(values[2]);
			probDeviation.push_back(values[3]);
		}
	}
	return handleInputData();	
}


// deal with the input data just read in from the file, if the file was read in correctly, return true.
// if return false, the input file has format errors
bool WeatherData::handleInputData()
{
	if(xCoors.empty())
		return false;
	// The total number of weather cells
	numPoints = xCoors.size();
	maxAlt = minAlt = (int)altitudes[0];
	maxProbDev = minProbDev = probDeviation[0];
	// get the range of each data components from the vectors
	bool inputFormatError = false;				// no error supposedly
	for(unsigned int i=0; i<numPoints; i++)
	{
		if(xCoors[i]>360 || xCoors[i]<-360 || yCoors[i]>360 || yCoors[i]<-360)
		{
			inputFormatError = true;
			break;
		}
#if defined(DEBUG)
		std::cout << xCoors[i] << " " << yCoors[i] << std::endl;
#endif
		// normalize the range into (-180, 180]
		if(xCoors[i]>180)
		{
			xCoors[i]-=360;
		}
		if(yCoors[i]>180)
		{
			yCoors[i]-=360;
		}
		if(altitudes[i]<minAlt)
		{
			minAlt = (int) altitudes[i];
		}
		if(altitudes[i]>maxAlt)	
		{
			maxAlt = (int) altitudes[i];
		}
		if(probDeviation[i]<minProbDev)	
		{
			minProbDev = probDeviation[i];
		}
		if(probDeviation[i]>maxProbDev)
		{
			maxProbDev = probDeviation[i];
		}
	}
	if(inputFormatError || minProbDev<0 || maxProbDev>1) 
	{
#if defined(DEBUG)
		std::cout << "minProbDev: " << minProbDev << "  maxProbDev: " << maxProbDev << std::endl;
#endif
		std::cerr << "\nWeather Data Content Error."<<std::endl;
		return false;									// the probability of deviation must be between 0 and 1, a messagebox is popped up
	}
	return true;										// the format is good and we are ready to draw
}

// a helper function used when reading in a file to help judge if the file format is correct
// test if we are going out of the end of the file, which means input format error. If format error, pop up a message box
bool WeatherData::testIndex(const int* readingIndex, const int *fileSize)
{
	if(*readingIndex < *fileSize)
	{
		return false;				// meaning there is not a problem
	}
	std::cerr << "\nFile Format Error!"<<std::endl;		// prompt that the file has format errors
	reset();
	return true;
}

// convert weather cell lati/long values to screen coordinates x/y values. 
// lati/long pair (cX, cY) is at x, y position (0, 0), latiPerPix, longPerPix are the lati/long units that a pixel represents
// after the conversion, compute the cellWidth/cellHeight value
void WeatherData::convertLatiLongHeightToXY(double cX, double cY, double latiPerPix, double longPerPix, double weatherCellWidth)
{
	for(unsigned int i=0; i<numPoints; i++)
	{
		
		xCoors[i] = (xCoors[i]-cX)/latiPerPix;
		yCoors[i] = (yCoors[i]-cY)/longPerPix;
    
		// altitude values usually vary from 10000 to 30000, difference is 5000
		// 0 is drawn on z=ALTITUDE_AT_BASE_PLANE; 15000 is at z=(15000-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL, etc.
		altitudes[i] = (altitudes[i]-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL;
	}
	cellHeight = (maxAlt - minAlt)/ALTITUDE_PER_PIXEL/4;// 10 pixels for 5000 feet altitude difference, suppose there are always 4 levels of weather data
	cellWidth = weatherCellWidth;									// a temporary value set for the width of each weather cell
}

// functions to return the deviation threshold values
double WeatherData::getMaxDevThres() const
{
	return maxProbDev;
}

double WeatherData::getMinDevThres() const
{
	return minProbDev;
}

// the number of weather cells in the weather data
unsigned int WeatherData::size() const
{
	return numPoints;
}

// set and get the probability value attached to the current weather data (the probability that the current weather happens)
void WeatherData::setProbability(double prob)
{
	if(prob>=0 && prob<=1)
	{
		probability = prob;
	}
	else if(prob < 0)
	{
		probability = 0;
	}
	else
	{
		probability = 1;
	}
}

double WeatherData::getProbability() const
{
	return probability;
}

// return the minimum altitude of the current weather data: the z coordinate in screen coordinate system
double WeatherData::getMinAlt() const
{
	return (minAlt-ALTITUDE_AT_BASE_PLANE) / ALTITUDE_PER_PIXEL;
}

double WeatherData::getMaxAlt() const
{
	return (maxAlt-ALTITUDE_AT_BASE_PLANE) / ALTITUDE_PER_PIXEL;
}

// given an index, return the corresponding weather cell's data, including its (x, y, z) and probability of deviation 
// as well as the cell's width and height(these 2 parameters are fixed among all cells)
bool WeatherData::getCellData(unsigned int index, double* x, double* y, double* alt, double* probDev, double* cWidth, double* cHeight) const
{
	if(index >= numPoints || index<0)
	{
		return false;
	}
	*x = xCoors[index];
	*y = yCoors[index];
	*alt = altitudes[index];
	*probDev = probDeviation[index];
	*cWidth = cellWidth;
	*cHeight = cellHeight;
	return true;
}
