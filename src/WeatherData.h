// Weather file has each line 4 elements: (x position, y position, altitude, probability of deviation)
// relative to a fixed point. We call each line a point and its related information
// the file ends with a '\n'

#pragma once

#ifndef WEATHERDATA_H
#define WEATHERDATA_H

#include <vector>
#include <string>

#define DRAW_STYLE_FILL 1
#define DRAW_STYLE_FRAME 2
#define DRAW_STYLE_3D 3

// the altitude at opengl coordinate system base plane z=0 is also altitude 0 feet from ground
#define ALTITUDE_AT_BASE_PLANE 5000
#define ALTITUDE_PER_PIXEL 500	// always use 10 pixels in z direction to denote 500 feet in the actual airspace

class WeatherData
{
public:
	WeatherData(void);
	~WeatherData(void);
	bool readInFileData(std::string fileName, double rangeMinLati, double rangeMinLong, double rangeMaxLati, double rangeMaxLong, double deviationThreshold);
	void reset();
  /**
    \brief Returns true iff there are no weather cells in this object
  */
	double getMaxDevThres() const;
	double getMinDevThres() const;
	unsigned int size() const;
	void setProbability(double prob);
	double getProbability() const;
	bool getCellData(unsigned int index, double* x, double* y, double* alt, double* probDeviation, double* cWidth, double* cHeight) const;
	void convertLatiLongHeightToXY(double cX, double cY, double latiPerPix, double longPerPix, double weatherCellWidth);
	double getMinAlt() const;
	double getMaxAlt() const;
private:
	bool handleInputData();			// in case the input file has format problem
	bool testIndex(const int* readingIndex, const int* fileSize);				// test if we are going out of the range of the file, meaning file format error
private:
	unsigned int numPoints;					// the number of points, hence weather cells(squares)
	std::vector<double> xCoors;			// the std::vector to store x positions of the weather cells
	std::vector<double> yCoors;			// the std::vector to store y positions of the weather cells
	std::vector<double> altitudes;		// the std::vector to store altitudes of the weather cells
	std::vector<double> probDeviation;	// the std::vector to store probability of deviation of the weather cells
	int minAlt;
	int maxAlt;				// the range of altitude
	double minProbDev;
	double maxProbDev;	// the range of deviation probabilities
	double cellWidth;				// when drawing, the width of each cell (the side length of a weather grid)
	double cellHeight;				// these 2 variables could be revised to std::vectors for each weather cell
private:
	double probability;				// the probability that the CURRENT weather data object is going to happen
public:
  static unsigned int readEnsembleMemberIndex( std::string fileName );
  static unsigned int readNumberOfEnsembleMembers( std::string fileName );
};

#endif
