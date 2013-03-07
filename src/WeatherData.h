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

using namespace std;

class WeatherData
{
public:
	WeatherData(void);
	~WeatherData(void);
	bool readInFileData(std::string fileName, double rangeMinLati, double rangeMinLong, double rangeMaxLati, double rangeMaxLong);
	void reset();
	double getMaxDevThres() const;
	double getMinDevThres() const;
	int size() const;
	void setProbability(double prob);
	double getProbability() const;
	bool getCellData(int index, float* x, float* y, float* alt, float* probDeviation, float* cWidth, float* cHeight) const;
	void convertLatiLongHeightToXY(double cX, double cY, double latiPerPix, double longPerPix);
	double getMinAlt() const;
	double getMaxAlt() const;
private:
	bool handleInputData();			// in case the input file has format problem
	bool testIndex(const int* readingIndex, const int* fileSize);				// test if we are going out of the range of the file, meaning file format error
private:
	int numPoints;					// the number of points, hence weather cells(squares)
	vector<float> xCoors;			// the vector to store x positions of the weather cells
	vector<float> yCoors;			// the vector to store y positions of the weather cells
	vector<float> altitudes;		// the vector to store altitudes of the weather cells
	vector<float> probDeviation;	// the vector to store probability of deviation of the weather cells
	int minAlt, maxAlt;				// the range of altitude
	float minProbDev, maxProbDev;	// the range of deviation probabilities
	double cellWidth;				// when drawing, the width of each cell (the side length of a weather grid)
	double cellHeight;				// these 2 variables could be revised to vectors for each weather cell
private:
	double probability;				// the probability that the CURRENT weather data object is going to happen
};

#endif
