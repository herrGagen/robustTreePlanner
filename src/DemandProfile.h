#ifndef DEMANDPROFILE_H
#define DEMANDPROFILE_H

#define PI 3.14159265

// for read in demand profiles, we always use 200 pixels to display the demand circle
#define TOTAL_DEMAND_CIRCLE_RADIUS_PIXEL 100
// the total number of entry nodes for each quadrant
#define NUM_ENTRY_NODES_PER_QUADRANT 8

#include <vector>
#include <string>

using namespace std;

class DemandProfile
{
public:
	DemandProfile(void);
	virtual ~DemandProfile(void);
	
public:
	bool readInFile(std::string fileName);
	void reset();
	void getDemandInfo(double* cX, double* cY, double* laPerPixel, double* loPerPixel, std::string* timeStart, std::string* timeEnd);
	void getRange(double* minLati, double* minLong, double* maxLati, double* maxLong);
	void generateDemandVector(std::vector<double> &demandRNPs, double startingAngle, double angularWidth, double nmilesPerPixel);
	
	// most x,y pairs we are talking about in this class are actually lati/long pairs
private:
	std::vector<double> xCoors;	/**< the x coordinates of the aircraft (latitudes) */
	std::vector<double> yCoors;	/**< the y coordinates of the aircraft (longitudes) */
	std::vector<double> rnp; /**< the rnp of each aircraft */
	std::string timeStart; /**< Simulation starting time */
        std::string timeEnd;   /**< Simulation ending time */
	double centerX;	// the x coordinate of the airport(latitude)
	double centerY;	// the y coordinate of the airport(longitude)
	double latiPerPixel; // conversion from latitude to pixel in x direction
	double longPerPixel; // conversion from longitude to pixel in y direction
	double altPerPixel;  
  double minX, minY, maxX, maxY;	// the ranges of the lati/long values
	unsigned int numDemands;		// the number of demands that are listed in the input file
private:
	bool handleInputData();	 // after reading in, deal with the input
        // test if we are going out of the range of the file, meaning file format error
	bool testIndex(const int* readingIndex, const int* fileSize);	
};

#endif
