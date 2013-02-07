#include "DemandProfile.h"
#include <fstream>
#include "math.h"
#include <iostream>
#include <cstring>
#include <cstdlib>

DemandProfile::DemandProfile(void)
{
	numDemands = 0;
	minX = maxX = minY = maxY = 0;
}

DemandProfile::~DemandProfile(void)
{
	xCoors.clear();
	yCoors.clear();
	rnp.clear();
}

// read in the demand profile from file, the content is passed in by the character array buffer, and the size is fileSize
bool DemandProfile::readInFile(char* buffer, int fileSize)			
{
	int wordStartingIndex = 0;											// record the start of a new word in the buffer
	int readingIndex = 0;												// record the position in the file, from 1 to fileSize											
	// the first step tries to read in the center location(airport location)
	while ( readingIndex < fileSize )
	{
		// the 4 separators: comma, line change('\n') and space. 13 "Carriage Return" is used together with '\n'
		while (buffer[readingIndex]!=',' && buffer[readingIndex]!='\n' && buffer[readingIndex]!=' ' && buffer[readingIndex]!=13)					
		{
			readingIndex++;
			if(testIndex(&readingIndex, &fileSize))						// if we are going out of range
				return false;
		}
		buffer[readingIndex] = '\0';									// the current char is a separator', mark it as the end of a string
		if( strcmp("SIMULATION_TIME_RANGE", &buffer[wordStartingIndex])==0)
		{
			wordStartingIndex = readingIndex+1;
			while(buffer[readingIndex]!=' ')
			{
				readingIndex++;
				if(testIndex(&readingIndex, &fileSize))
					return false;
			}
			buffer[readingIndex] = '\0';
			string tempStartingTime(&buffer[wordStartingIndex]);		// the current word is the starting time of the demand time range
			timeStart = tempStartingTime;								// store it in the timeStart string
			// next word would be the ending time of the demand time range
			wordStartingIndex = readingIndex+1;
			while(buffer[readingIndex]!='\n')
			{
				readingIndex++;
				if(testIndex(&readingIndex, &fileSize))
					return false;
			}
			buffer[readingIndex] = '\0';
			string tempEndingTime(&buffer[wordStartingIndex]);			// the current word is the starting time of the demand time range
			timeEnd = tempEndingTime;									// store it in the timeStart string
		}
		if( strcmp("CENTER_LOCATION", &buffer[wordStartingIndex])==0)	// if we find the center location "CENTER_LOCATION"
		{
			// the center location is revealed in the next strings
			wordStartingIndex = readingIndex+1;							// read in the x position of the center first
			while(buffer[readingIndex]!=' ')
			{
				readingIndex++;
				if(testIndex(&readingIndex, &fileSize))					// if we are going out of range
					return false;
			}
			buffer[readingIndex] = '\0';								// mark it as the end of a word
			centerX = atof(&buffer[wordStartingIndex]);
			wordStartingIndex = readingIndex+1;							// next read in the y position of the center
			while(buffer[readingIndex]!='\n')
			{
				readingIndex++;
				if(testIndex(&readingIndex, &fileSize))					// if we are going out of range
					return false;
			}
			buffer[readingIndex] = '\0';								// mark it as the end of a word
			centerY = atof(&buffer[wordStartingIndex]);
			wordStartingIndex = readingIndex+1;
			break;
		}
		wordStartingIndex = readingIndex+1;								// read the next word
	}
	// the second step tries to find the start position of the flight information
	while (readingIndex < fileSize)
	{
		while(buffer[readingIndex]!=',' && buffer[readingIndex]!='\n' && buffer[readingIndex]!=' ' && buffer[readingIndex]!=13)
		{
			readingIndex++;
			if(testIndex(&readingIndex, &fileSize))						// if we are going out of range
				return false;
		}
		buffer[readingIndex] = '\0';
		// in between "BEGIN_FLIGHTS" and "END_FLIGHTS" are the flight information
		if( strcmp("BEGIN_FLIGHTS", &buffer[wordStartingIndex])==0)		// the current word is "BEGIN_FLIGHTS"
		{
			wordStartingIndex = readingIndex+1;							// set the starting position of the next word
			break;
		}
		wordStartingIndex = readingIndex+1;
	}
	// the third step reads in the flight arrival infomation
	while(readingIndex<fileSize)
	{
		for(int i=0; i<5; i++)											// the first 5 elements of each line are irrelavant
		{
			while(buffer[readingIndex]!=',')							// find the first 5 commas
			{
				readingIndex++;
				if(testIndex(&readingIndex, &fileSize))					// if we are going out of range
					return false;
			}
			readingIndex++;
		}
		// after the first 5 elements, it comes the rnp information that we are going to store in vector<float> rnp
		wordStartingIndex = readingIndex;
		while(buffer[readingIndex]!=',')
		{
			readingIndex++;
			if(testIndex(&readingIndex, &fileSize))
				return false;
		}
		buffer[readingIndex] = '\0';
		rnp.push_back(atof(&buffer[wordStartingIndex]));
		wordStartingIndex = ++readingIndex;								// reading index is at the start of the next word
		// the 6th to the 9th elements are irrelevant, we just skip these 3 elements by skipping 3 commas
		for(int i=0; i<3; i++)
		{
			while(buffer[readingIndex]!=',')
			{
				readingIndex++;
				if(testIndex(&readingIndex, &fileSize))
					return false;
			}
			readingIndex++;
		}
		// after the first 9 elements, it comes the information that we need: the flight (x, y) coordinates
		wordStartingIndex = readingIndex;								// the start of the next word
		while(buffer[readingIndex]!=',')								// read in the x coordinate in to xCoors
		{
			readingIndex++;
			if(testIndex(&readingIndex, &fileSize))						// if we are going out of range
				return false;
		}
		buffer[readingIndex] = '\0';
		xCoors.push_back(atof(&buffer[wordStartingIndex]));				// store in xCoors
		wordStartingIndex = readingIndex+1;
		while(buffer[readingIndex]!='\n')								// read in the y coordinate in to yCoors(separated by line change'\n')
		{
			readingIndex++;
			if(testIndex(&readingIndex, &fileSize))						// if we are going out of range
				return false;
		}
		buffer[readingIndex] = '\0';
		yCoors.push_back(atof(&buffer[wordStartingIndex]));				// store in yCoors
		wordStartingIndex = readingIndex+1;
		if(buffer[wordStartingIndex] == 'E')							// if the next word is "END_FLIGHTS", then we are done reading
			break;
	}
	return handleInputData();
}

bool DemandProfile::handleInputData()
{
	if(centerX>360 || centerY>360)							// the lati/long of the center is invalid
		return false;
	if(centerX>180)	centerX -=360;							// normalize the lati/longs to (-180, 180] range
	if(centerY>180)	centerY -=360;
	numDemands = xCoors.size();								// the total number of demands(flights) in the profile
	// the demand profile is almost on a circle centers at position(centerX, centerY)
	minX = maxX = centerX;
	minY = maxY = centerY;
	// check if there is input content error, and regulate the lati/long to the range of (-180, 180]
	for(int i=0; i<numDemands; i++) 
	{
		if(rnp[i]<0 || xCoors[i]>360 || yCoors[i]>360)
		{
			cerr<<"\nFile Content Error!"<<endl;				// prompt that the file has format errors
			return false;
		}
		if(xCoors[i]>180)									// change the format to (-180, 180]
			xCoors[i]-=360;
		if(yCoors[i]>180)
			yCoors[i]-=360;
		if(xCoors[i]<minX)	minX = xCoors[i];
		if(yCoors[i]<minY)	minY = yCoors[i];
		if(xCoors[i]>maxX)	maxX = xCoors[i];
		if(yCoors[i]>maxY)	maxY = yCoors[i];
	}
	// find the latitude diffference
	// the demand circle has total diameter of 200 pixels(a given value, best value chosen for display purpose) 
	latiPerPixel = (maxX - minX)/(2*TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL);						
	longPerPixel = (maxY - minY)/(2*TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL);
	return true;
}

void DemandProfile::reset()
{
	xCoors.clear();
	yCoors.clear();
	rnp.clear();
	numDemands = 0;
	minX = maxX = minY = maxY = 0;
}

// test if we are going out of the end of the file, which means input format error. Pop up a messagebox to prompt format error and reset the demandprofile
bool DemandProfile::testIndex(const int* readingIndex, const int *fileSize)
{
	if(*readingIndex < *fileSize)
		return false;						// meaning there is not a problem
	cerr<<"\n File Format Error!"<<endl;	// prompt that the file has format errors
	reset();
	return true;							// means we are out of file total size
}

// thsi function is called when reading in weather data. The function is used to do data conversion, from
// the read in weather cell in lati/long to their x/y positions. 
void DemandProfile::getDemandInfo(double* cX, double* cY, double* laPerPixel, double* loPerPixel, string* timeS, string* timeE)
{
	*cX = centerX;
	*cY = centerY;
	*laPerPixel = latiPerPixel;
	*loPerPixel = longPerPixel;
	*timeS = timeStart;
	*timeE = timeEnd;
}

// The function is also used to do weather data triming,
// where the data cells that are way too far from the demand are not stored in the memory at all
void DemandProfile::getRange(double *minLati, double *minLong, double *maxLati, double *maxLong)
{
	*minLati = minX;	*minLong = minY;
	*maxLati = maxX;	*maxLong = maxY;
}

// generate a set of 8 entry points and their rnp values based on the starting angle of the quadrant
void DemandProfile::generateDemandVector(std::vector<float> &demand, double startingAngle, double endingAngle, float nmilesPerPixel)
{
	float demandRNPS[NUM_ENTRY_NODES_PER_QUADRANT] = {0};
	double angleIncrement = (endingAngle - startingAngle)/NUM_ENTRY_NODES_PER_QUADRANT;
	for(int i=0; i<numDemands; i++)
	{
		// compute the x, y coordinates in the opengl coordinate system for each demand, where the center is (0, 0)
		double tempX = (xCoors[i]-centerX)/latiPerPixel;	
		double tempY = (yCoors[i]-centerY)/longPerPixel;
		if(tempX!=0)
		{
			double angle = atan(tempY/tempX);					// the center angle of the point
			if(tempX<0)											// if the point is in the 2nd or the 3rd quadrant
				angle = angle+PI;
			else if(tempX>0 && tempY<0)							// in the 4th quadrant
				angle = angle+2*PI;								// make it in the 0--2*PI range
			// now we have the angle with in range (0, 2*PI)
			if(endingAngle<2*PI)								// if the quadrant's angle
			{
				if(angle>startingAngle && angle<=endingAngle)
				{
					int tempIndex = int((angle-startingAngle)/angleIncrement);
					demandRNPS[tempIndex] = max(demandRNPS[tempIndex], rnp[i]);	 // compute the max rnp value for this range of entry node
				}
			}
			else												// starting and ending angle are on different sides of the positive x axis
			{
				if(angle>startingAngle || angle<endingAngle-2*PI)
				{
					int tempIndex = angle>=startingAngle ? int((angle-startingAngle)/angleIncrement) :  int((angle+2*PI-startingAngle)/angleIncrement);
					demandRNPS[tempIndex] = max(demandRNPS[tempIndex], rnp[i]);	 // compute the max rnp value for this range of entry node
				}
			}			
		}
	}
	// convert to opengl coordinate system and in the project, rnp is half the lane width
	for(int i=0; i<NUM_ENTRY_NODES_PER_QUADRANT; i++)		
		demand.push_back(demandRNPS[i]/(2*nmilesPerPixel));						// store them into the resulting vector
}
