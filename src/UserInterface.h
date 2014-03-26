#include "Quadrant.h"
#include "WeatherData.h"
#include "DemandProfile.h"
#include "RoutingDAG.h"
#include <string>
#include "InputFileReader.h"

#if defined(SUPPRESS_OUTPUT)
#define cout ostream(0).flush()
#endif

#define QUADRANT_NOT_GENERATED 0			// define if the quadrant is generated or not
#define QUADRANT_GENERATED 1

#define WEATHER_NOT_READ_IN 0				// define if the weather data files are read in
#define WEATHER_READ_IN 1

#define DEMAND_NOT_READ_IN 0				// define if the demand profile file is read in
#define DEMAND_READ_IN 1

#define ROUTINGDAG_NOT_GENERATED 0			// define if the routing DAG is generated
#define ROUTINGDAG_GENERATED 1

#define OPER_FLEX_GENERATED 0				// define if the operational flexibility is generated
#define OPER_FLEX_NOT_GENERATED 1

#define PI 3.14159265						// the value of PI

#ifndef USERINTERFACE_H
#define USERINTERFACE_H

class UserInterface
{
public:
	UserInterface();						// constructor & destructor
	~UserInterface();
	/****************************************************************************************************************************/
	// general purpose functions
public:
	void reset();							// reset the status to a brand new routing instance
	bool ProgramBegins(std::string inputFile);					// the project starts excuting from this function
	void saveTreeInformation();				// after generating the tree, save the information into an .xml file
	void inputOperationalFlexibility();		// ask user to input parameters related to operational flexity
private:
	void resetHelper();
	/****************************************************************************************************************************/
	// get user/file input or save to files functions
 private:
	bool makeRTPTreeAndFinish();
	bool makeDAG();
	bool editQuadrant();					// ask for user input and edit the quadrant information
	bool readWeatherData();					// read in from weather data files 
	bool readDemandProfile();				// read in the demand profile information
	bool inputDemand();						// ask the user to input demand information
	void printQuadrantAndDemandInfo();		// use standard output to print the quadrant information on the screen
	bool inputDemandValid(std::string &input, int numDemands);		// tell if the user input demand is valid or not
	/****************************************************************************************************************************/
	// tree generating related functions
	bool generateTree();					// generate a bottommost robust tree
	bool tautenTree();						// tauten the generated bottommost tree
	/****************************************************************************************************************************/
private:
	InputFileReader inputs;
	std::string outputFileName;
        // the quadrant in the routing context
	Quadrant* quadrant;				
	// a set of weather data files
	std::vector<WeatherData> weatherDataSets;
        // demand profile information
	DemandProfile* demandProfile;
        // key component, the routing DAG structure
	RoutingDAG* routingDAG;	
        // the std::vector storing RNP values for each entry node
	std::vector<double> demandRNPs;
	// the threshold that a weather cell is considered hazardous
	double deviationThreshold;				
	// the threshold that we consider a thick edge or disk to be clear of weather
	double nodeEdgeThreshold;				
public:
	const std::vector<double> &getDemandRNPs() const{ return demandRNPs; }
	double getDeviationThreshold() const{ return deviationThreshold; }
	double getNodeEdgeThreshold() const { return nodeEdgeThreshold; }
	double getCenterLati() const { return centerLati; }
	double getCenterLong() const { return centerLong; }
	double getLatiPerPixel() const { return latiPerPixel; }
	double getLongPerPixel() const { return longPerPixel; }
	RoutingDAG *getRoutingDAG() const { return routingDAG; }
	const std::vector<WeatherData> &getWeatherDataSets() const {return weatherDataSets;}
	std::string getOutputFileName() const { return outputFileName; }
	void setOutputFileName( std::string newName) { outputFileName = newName; }
    const std::vector<double> &getOperFlex() const { return inputs.getOperFlex(); }
	/********************************/
	// control variables
private:
        // mark if a quadrant was generated already
	int ctrl_QuadGenerated;
        // mark if a weather data file is already read in
	int ctrl_WeatherReadIn;		
        // if the demand was read in already
	int ctrl_DemandReadIn;
        // if a routingDAG is generated
	int ctrl_RoutingDAGGenerated;	
        // if the operational flxibility pairs are generated for an edge or a node
	int ctrl_OperFlexGenerated;	
	/******************************/
	// variables related to coordinate system conversion between 
        // OPENGL(used in the windows version for displaying) and lat/lon
private: 
        // the lati/long value of the center point (0, 0)
	double centerLati,centerLong;	
	// lati/long that a single pixel represents
	double latiPerPixel, longPerPixel;	
        // the time range of simulation
	std::string startTime;
        std::string endTime;		

  /******* Functions to help us only use relevant weather points *************/
  static std::pair<double, double> latLonOfPointAlongLineWithBearing(double lat, double lon, double dist, double bearing);
  void calculateBoundingBox(double &minLat, double &minLon, double &maxLat, double &maxLon) const;

};

#endif

