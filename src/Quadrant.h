#ifndef QUANDRANT_H
#define QUANDRANT_H

#include "WeatherData.h"
#include "RoutingDAG.h"

#define PI 3.14159265

#define QUADRANT_PLANE 1				// constants describing if there is height given to the quadrant
#define QUADRANT_LIFTED 2

class Quadrant
{
public:
	Quadrant(double cX = 0, double cY = 0, double ang = PI/6, double iR = 10, double oR = 35, double iH = 0, double oH = 0);
	virtual ~Quadrant(void);
	void setcX(double cX);
	void setcY(double cY);
	virtual void setAngle(double ang);
	void setiRadius(double iR);
	void setoRadius(double oR);
	void setiHeight(double iH);
	void setoHeight(double oH);
	double getcX();
	double getcY();
	double getAngle();
	double getiRadius();
	double getoRadius();
	double getiHeight();
	double getoHeight();
	int getLiftStatus();
	void setLiftStatus(int status);
	void setQuadrant(double cX = 0, double cY = 0, double ang = PI/6, double iR = 10, double oR = 35, double iH = 0, double oH = 0);
	bool demandFeasible(std::vector<double> &rnps);
	void reset();
	/*********************************************************************************************/
	// functions used to generate routing graph structures
	bool generateDAG(std::vector<double> rnps, int n, double effectiveThres, double routingThres, const std::vector<WeatherData> &wData, RoutingDAG* rDAG, double qAngleOffset, int numFixNodes);
	/*********************************************************************************************/
	
private:
	int liftStatus;		// if the quadrant is lifted on not
	double centerX;		// the x coordinate of the center
	double centerY;		// the y coordinate of the center
	double angle;		// the angle relative to the x axis, from 0 to 2*PI
	double iRadius;		// inner radius of the quadrant
	double oRadius;		// outer radius of the quadrant
	double liftediRadius;	// the inner and outer radius after the quadrant is lifted
	double liftedoRadius;
	double iHeight;		// the inner height of the boundary of the quadrant
	double oHeight;		// the outer height of the boundary of the quadrant
	double cHeight;		// the height of the center of the quadrant
	
private:
	/*********************************************************************************************/
	// functions used to generate routing graph structures
	bool generateEntryAndFixNodes(std::vector<double> rnps, double effectiveThres, double routingThres, const std::vector<WeatherData> &wData, RoutingDAG* rDAG, double quadAngleOffset, unsigned int numFixNodes);
	void generateRoutingDAGInternalNodes(RoutingDAG* rDAG, std::vector<double> rnps, double quadAngleOffset);
	/*********************************************************************************************/
};

#endif
