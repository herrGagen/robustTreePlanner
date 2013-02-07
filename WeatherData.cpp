#include "WeatherData.h"
#include "math.h"
#include <algorithm>
#include "RoutingDAG.h"
#include <time.h>
#include <iostream>
#include <cstring>
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
bool WeatherData::readInFileData(char* buffer, int fileSize, double rangeMinLati, double rangeMinLong, double rangeMaxLati, double rangeMaxLong)
{
	/******************************************************************************************************/
	// test file read in time
	/*LARGE_INTEGER m_lpFrequency, lpPerformanceCount2,lpPerformanceCount1;  //variables for time counting
	 QueryPerformanceFrequency(&m_lpFrequency);
	 QueryPerformanceCounter(&lpPerformanceCount1);*/
	/******************************************************************************************************/
	reset();											// first reset the vectors to empty, then read in
	// Starting here, we deal with the input data
	int wordStartingIndex = 0;							// record the start of a new word in the buffer
	int readingIndex = 0;								// record the position in the file, from 1 to fileSize
	while ( readingIndex < fileSize )
	{
		while (buffer[readingIndex]!=' ' && buffer[readingIndex]!='\n')
			readingIndex++;
		buffer[readingIndex]='\0';						//mark it as the end of a string
		// read out the probability of the current weather file
		if(strcmp("Probability", &buffer[wordStartingIndex])==0)
		{
			wordStartingIndex = readingIndex+1;
			while(buffer[readingIndex]!=' ' && buffer[readingIndex]!='\n')
				readingIndex++;
			buffer[readingIndex]='\0';
			probability = atof(&buffer[wordStartingIndex]);
			wordStartingIndex = readingIndex+1;
			break;
		}
		wordStartingIndex = readingIndex+1;
	}
	// after reading the probability, read in the data file content
	float tempX, tempY;									// read in these 2 values first to see if they are with in the range of x and y	
	bool pushIn = true;									// denote if the current 4-pair is valid based on the range
	while(readingIndex < fileSize)
	{
		while (buffer[readingIndex]!=',')				// the first 3 elements of each line are separated by comma
			readingIndex++;
		buffer[readingIndex] = '\0';					// the current char is ',', mark it as the end of a string
		tempX = atof(&buffer[wordStartingIndex]);
		wordStartingIndex = readingIndex+1;				// set to the start of the next C style string
		// do the same thing for the next input number, write the codes 3 times to save time
		while (buffer[readingIndex]!=',')	// the first 3 elements of each line are separated by comma
			readingIndex++;
		buffer[readingIndex] = '\0';					// the current char is ',', mark it as the end of a string
		tempY = atof(&buffer[wordStartingIndex]);
		/************************************************************************************************/
		// now we have the (tempX, tempY) pair and we do the range query to see if this is valid a pair
		// if they are with in the range, save the information into the corresponding vectors
		pushIn = false;
		if(tempX>rangeMinLati && tempX<rangeMaxLati && tempY>rangeMinLong && tempY<rangeMaxLong)
		{
			xCoors.push_back(tempX);	
			yCoors.push_back(tempY);
			pushIn = true;
		}
		/************************************************************************************************/
		wordStartingIndex = readingIndex+1;				// set to the start of the next C style string
		// do the same thing for the next input number, write the codes 3 times to save time
		while (buffer[readingIndex]!=',')	// the first 3 elements of each line are separated by comma
			readingIndex++;
		buffer[readingIndex] = '\0';					// the current char is ',', mark it as the end of a string
		if(pushIn)										// if the lati/long pair is valid
			altitudes.push_back(atof(&buffer[wordStartingIndex]));
		wordStartingIndex = readingIndex+1;				// set to the start of the next C style string
		// for the last element, the probability of deviation, its separated by '\n'
		// do the same thing for the next input number, write the codes 3 times to save time
		while (buffer[readingIndex]!='\n')	// the last element of each line is separated by '\n', so is the last element of the file
			readingIndex++;
		buffer[readingIndex] = '\0';					// the current char is ',', mark it as the end of a string
		if(pushIn)										// if the lati/long pair is valid
			probDeviation.push_back(atof(&buffer[wordStartingIndex]));
		wordStartingIndex = readingIndex+1;				// set to the start of the next(possible) C style string
		readingIndex++;									// this is only for the last line, where we are jump out of the while
	}
	/******************************************************************************************************/
	// end of testing file read in time
	/*QueryPerformanceCounter(&lpPerformanceCount2);                                // time query
	 double time = (double)(lpPerformanceCount2.LowPart - lpPerformanceCount1.LowPart)/m_lpFrequency.LowPart;
	 TRACE("The Time Needed is %g\n", time);*/
	/******************************************************************************************************/
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
	maxAlt = minAlt = altitudes[0];
	maxProbDev = minProbDev = probDeviation[0];
	// get the range of each data components from the vectors
	bool inputFormatError = false;				// no error supposedly
	for(int i=0; i<numPoints; i++)
	{
		if(xCoors[i]>360 || xCoors[i]<-360 || yCoors[i]>360 || yCoors[i]<-360)
		{
			inputFormatError = true;
			break;
		}
		// normalize the range into (-180, 180]
		if(xCoors[i]>180)	xCoors[i]-=360;
		if(yCoors[i]>180)	yCoors[i]-=360;
		if(altitudes[i]<minAlt)	minAlt = altitudes[i];
		if(altitudes[i]>maxAlt)	maxAlt = altitudes[i];
		if(probDeviation[i]<minProbDev)	minProbDev = probDeviation[i];
		if(probDeviation[i]>maxProbDev)	maxProbDev = probDeviation[i];
	}
	if(inputFormatError || minProbDev<0 || maxProbDev>1) 
	{
		cerr<<"\nWeather Data Content Error."<<endl;
		return false;									// the probability of deviation must be between 0 and 1, a messagebox is popped up
	}
	return true;										// the format is good and we are ready to draw
}

// a helper function used when reading in a file to help judge if the file format is correct
// test if we are going out of the end of the file, which means input format error. If format error, pop up a message box
bool WeatherData::testIndex(const int* readingIndex, const int *fileSize)
{
	if(*readingIndex < *fileSize)
		return false;						// meaning there is not a problem
	cerr<<"\nFile Format Error!"<<endl;		// prompt that the file has format errors
	reset();
	return true;
}

// convert weather cell lati/long values to screen coordinates x/y values. 
// lati/long pair (cX, cY) is at x, y position (0, 0), latiPerPix, longPerPix are the lati/long units that a pixel represents
// after the conversion, compute the cellWidth/cellHeight value
void WeatherData::convertLatiLongHeightToXY(double cX, double cY, double latiPerPix, double longPerPix)
{
	for(int i=0; i<numPoints; i++)
	{
		xCoors[i] = (xCoors[i]-cX)/latiPerPix;
		yCoors[i] = (yCoors[i]-cY)/longPerPix;
		// altitude values usually vary from 10000 to 30000, difference is 5000
		// 0 is drawn on z=ALTITUDE_AT_BASE_PLANE; 15000 is at z=(15000-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL, etc.
		altitudes[i] = (altitudes[i]-ALTITUDE_AT_BASE_PLANE)/ALTITUDE_PER_PIXEL;
	}
	cellHeight = (maxAlt - minAlt)/ALTITUDE_PER_PIXEL/4;// 10 pixels for 5000 feet altitude difference, suppose there are always 4 levels of weather data
	cellWidth = 1.25;									// a temporary value set for the width of each weather cell
}

// functions to return the deviation threshold values
double WeatherData::getMaxDevThres()
{
	return maxProbDev;
}

double WeatherData::getMinDevThres()
{
	return minProbDev;
}

// the number of weather cells in the weather data
int WeatherData::size()
{
	return numPoints;
}

// set and get the probability value attached to the current weather data (the probability that the current weather happens)
void WeatherData::setProbability(double prob)
{
	if(prob>=0 && prob<=1)
		probability = prob;
}

double WeatherData::getProbability()
{
	return probability;
}

// return the minimum altitude of the current weather data: the z coordinate in screen coordinate system
double WeatherData::getMinAlt()
{
	return (minAlt-ALTITUDE_AT_BASE_PLANE) / ALTITUDE_PER_PIXEL;
}

double WeatherData::getMaxAlt()
{
	return (maxAlt-ALTITUDE_AT_BASE_PLANE) / ALTITUDE_PER_PIXEL;
}

// given an index, return the corresponding weather cell's data, including its (x, y, z) and probability of deviation 
// as well as the cell's width and height(these 2 parameters are fixed among all cells)
bool WeatherData::getCellData(int index, float* x, float* y, float* alt, float* probDev, float* cWidth, float* cHeight)
{
	if(index>=numPoints || index<0)
		return false;
	*x = xCoors[index];
	*y = yCoors[index];
	*alt = altitudes[index];
	*probDev = probDeviation[index];
	*cWidth = cellWidth;
	*cHeight = cellHeight;
	return true;
}
