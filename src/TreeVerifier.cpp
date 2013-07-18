#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

#include "RoutingDAG.h"
#include "TreeVerifier.h"

namespace
{
	/**
	\brief Returns the creat circle distance between two points on earth

	\parm lat1 Latitude of the first point
	\parm lon1 Longitude of the first point
	\parm lat2 Latitude of the second point
	\parm lon2 Longitude of the second point

	*/	
	double haversineDistance(double lat1, double lon1, double lat2, double lon2)
	{
		const double rEarthNm = 3443.89849; // Earth's radius in nmi
		double deltaLat = lat2-lat1; // difference in latitude
		double deltaLon = lon2-lon1; // difference in longitude

		double haversineDeltaLat = sin(deltaLat/2)*sin(deltaLat/2);
		double haversineDeltaLon = sin(deltaLon/2)*sin(deltaLon/2);
		double a = haversineDeltaLat + cos(lat1)*cos(lat2)*haversineDeltaLon;

		return rEarthNm * 2 * atan2(sqrt(a), sqrt(1-a));
	}

	/**
		\brief Finds the derivative of the haversine function at (lat1 lon1, lat2 lon2)

		Please note order of arguments.  

		This takes a numeical derivative of the haversine great circle distane at points
		[lat1, lon1] and [lat2, lon2] 
		Where [lat1 lon1] lie along the parameterized line (x0,y0) + t(dLat1, dLat2)
	*/
	double haversineDeriv(double lat1, double lon1, double dLat1, double dLon1, double lat2, double lon2)
	{
		const double epsilon = .0001;
		double dNegDelta = haversineDistance(lat1-dLat1*epsilon, lon1-dLon1*epsilon, lat2, lon2);
		double dDelta    = haversineDistance(lat1+dLat1*epsilon, lon1+dLon1*epsilon, lat2, lon2);

		return (dDelta - dNegDelta)/(2*epsilon);
	}
}

void TreeVerifier::appendReportToFiles(std::string invalidName, std::string validName)
{
	std::cout << "Verifying routing DAG is clear of weather." << std::endl;
	int badNodes = countInvalidNodes();
	int badEdges = countInvalidEdges();
	if(badNodes + badEdges > 0)
	{
		ofstream invalidOut;
		invalidOut.open(invalidName.c_str(),ios::app);
		invalidOut << ui.getOutputFileName() << " had ";
		invalidOut << badNodes << " bad nodes and " << badEdges << "bad edges." << std::endl;
		invalidOut.close();
		std::cout << "RoutingDAG not clear." << std::endl;
	}
	else
	{
		ofstream validOut;
		validOut.open(validName.c_str(), ios::app);
		validOut << ui.getOutputFileName() << std::endl;
		validOut.close();
		std::cout << "RoutingDAG clear of weather." << std::endl;
	}

}

TreeVerifier::TreeVerifier(const UserInterface &UI):
	ui(UI)
{
	centerLat = ui.getCenterLati();
	centerLon= ui.getCenterLong();
	latPerPixel = ui.getLatiPerPixel();
	lonPerPixel = ui.getLongPerPixel();
}
/**
\brief Counts number of nodes likely to be impacted by dangerous weather

Safe here is defined (per node) as:
For each weather set
If (any point in weather set is closer than laneWidth+cellWidth
AND that weather's deviationProbability > deviationThreshold)
Node's probability of being unsafe += weatherSet's probability

A node is safe if its probability of being safe > nodeEdgeThreshold

*/
unsigned int TreeVerifier::countInvalidNodes() const
{	

	unsigned int numInvalidNodes = 0;
	const RoutingDAG &dag = *(ui.getRoutingDAG() );
	for(unsigned int i = 0; i< dag.getNumNodes(); i++)
	{
		Node *thisNode = dag.getNodePointer(i);
		// If this node isn't in the tree, who cares if it impacts weather?
		if(!thisNode->isTreeNode() )
		{
			continue;
		}
		double pThisNodeIsSafe = totalProbNodeIsSafe(thisNode);
		if(pThisNodeIsSafe < ui.getNodeEdgeThreshold() )
		{
			numInvalidNodes++;
		}
	}
	return numInvalidNodes;
}

/**
	\brief Adds probabilit of existence for all weather ensembles this node avoids

	\param node The node whose safety we're verifying.
*/
double TreeVerifier::totalProbNodeIsSafe( Node *node ) const
{
	double pThisNodeIsSafe = 0;
	const std::vector<WeatherData> &weatherDataSets = ui.getWeatherDataSets();
	for( unsigned int j = 0; j<weatherDataSets.size(); j++)
	{
		if( doesNodeAvoidDangerousWeatherInEnsemble( node, j ) )
		{
			pThisNodeIsSafe += weatherDataSets[j].getProbability();
		}
	}
	return pThisNodeIsSafe;
}


/**
	\brief Checks is this node avoids dangerous weather in weather ensemble ensembleIndex

	Weather is dangerous if its deviationProbability is above ui.getDeviationThreshold()

	Weather is avoided if it's farther away than the weather's cell size + the RNP associated
	with this particular node

	\param Node node being checked
	\parm ensembleIndex Which ensemble in weatherDataSets we're checking for impacting our weather.

	\retval True is this node does not impact on dangerous weather in this ensemble
*/
bool TreeVerifier::doesNodeAvoidDangerousWeatherInEnsemble( Node *node, unsigned int ensembleIndex) const
{
	double deviationThresh = ui.getDeviationThreshold();
	const std::vector<WeatherData> &weatherDataSets = ui.getWeatherDataSets();
	const WeatherData &wData = weatherDataSets[ensembleIndex];

	// Storage for weather data's location
	// Don't shoot the messenger.
	double x;
	double y;
	double z;
	double cellWidth;
	double cellHeight;
	double deviationProbability;

	for(unsigned int i=0; i<wData.size() ; i++)
	{
		wData.getCellData(i, &x, &y, &z, &deviationProbability, &cellWidth, &cellHeight); 
		// If this weather isn't dangerous, don't bother finding its position.
		if(deviationProbability < ui.getDeviationThreshold() )
		{
			continue;
		}
		double safeDistance = node->getDrawingRNP() + cellWidth;
		// wLat/wLon weather's position
		double wLat = xToLat(x);
		double wLon = yToLon(y);
		// nLat / nLon node's position
		double nLat = xToLat( node->getX() );
		double nLon = yToLon( node->getY() );
		double thisDistance = haversineDistance(nLat, nLon, wLat, wLon );
		if(thisDistance < safeDistance)
		{
			return false;
		}
	}
	return true; // We only reach this if there was no impacting dangerous weather.

}

/**
	\brief Counts number of edges likely to be impacted by dangerous weather

	Safe here is defined (per edge) as:
	For each weather set
	If (any point in weather set is closer than laneWidth+cellWidth
	AND that weather's deviationProbability > deviationThreshold)
	Edge's probability of being unsafe += weatherSet's probability

	An edge is safe if its probability of being unsafe < 1-nodeEdgeThreshold

*/
unsigned int TreeVerifier::countInvalidEdges() const
{	

	unsigned int numInvalidEdges = 0;
	const RoutingDAG &dag = *(ui.getRoutingDAG() );
	for(unsigned int i = 0; i< dag.getNumEdges(); i++)
	{
		Edge *thisEdge = dag.getEdgePointer(i);
		// If this edge isn't in the tree, who cares if it impacts weather?
		if(!thisEdge->isTreeEdge() )
		{
			continue;
		}
		double pThisEdgeIsSafe = totalProbEdgeIsSafe(thisEdge);
		if(pThisEdgeIsSafe < ui.getNodeEdgeThreshold() )
		{
			numInvalidEdges++;
		}
	}
	return numInvalidEdges;
}

/**
	\brief Adds probability of existence for all weather ensembles this edge avoids

	\param edge The edge whose safety we're verifying.
*/
double TreeVerifier::totalProbEdgeIsSafe( Edge *edge ) const
{
	double pThisEdgeIsSafe = 0;
	const std::vector<WeatherData> &weatherDataSets = ui.getWeatherDataSets();
	for( unsigned int j = 0; j<weatherDataSets.size(); j++)
	{
		if( doesEdgeAvoidDangerousWeatherInEnsemble( edge, j ) )
		{
			pThisEdgeIsSafe += weatherDataSets[j].getProbability();
		}
	}
	return pThisEdgeIsSafe;
}


/**
	\brief Checks is this Edge avoids dangerous weather in weather ensemble ensembleIndex

	Weather is dangerous if its deviationProbability is above ui.getDeviationThreshold()

	Weather is avoided if cell size + the RNP associated with this particular Edge is greater
	than the distance to closest approach of the weather.

	\param Edge Edge being checked
	\parm ensembleIndex Which ensemble in weatherDataSets we're checking for impacting our weather.
*/
bool TreeVerifier::doesEdgeAvoidDangerousWeatherInEnsemble( Edge *edge, unsigned int ensembleIndex) const
{
	double deviationThresh = ui.getDeviationThreshold();
	const std::vector<WeatherData> &weatherDataSets = ui.getWeatherDataSets();
	const WeatherData &wData = weatherDataSets[ensembleIndex];

	// Storage for weather data's location
	// Don't shoot the messenger.
	double x;
	double y;
	double z;
	double cellWidth;
	double cellHeight;
	double deviationProbability;

	for(unsigned int i=0; i<wData.size() ; i++)
	{
		wData.getCellData(i, &x, &y, &z, &deviationProbability, &cellWidth, &cellHeight); 
		// If this weather isn't dangerous, don't bother finding its position.
		if(deviationProbability < ui.getDeviationThreshold() )
		{
			continue;
		}
		double safeDistance = edge->getDrawingRNP() + cellWidth;
		// wLat/wLon weather's position
		double wLat = xToLat(x);
		double wLon = yToLon(y);
		// nLat / nLon Edge's position

		if(boundDistWeatherToEdge(edge, wLat, wLon) < safeDistance)
		{
			double thisDistance = shortestDistWeatherToEdge(edge, wLat, wLon);
			if(thisDistance < safeDistance)
			{
				return false;
			}
		}
	}
	return true; // We only reach this if there was no impacting dangerous weather.

}

/**
	\brief Finds the closest that point (weatherLat, weatherLon) is to an edge

	because this uses great circle distance, this is non-trivial

	We are going to parameterize the edge as:
	(lat1, lon1) + t(lat2, lon2) where t \in [0,1]

	Then we look for the zeros of the partial derivative of distance w.r.t. t
	
	But we have no analytical form of that partial derivative, so we're going to 
	approximate it numerically.

*/
double TreeVerifier::shortestDistWeatherToEdge( Edge *edge, double weatherLat, double weatherLon ) const
{
	const double EPSILON = .01;

	double lat1 = xToLat( edge->getX1() );	
	double lon1 = yToLon( edge->getY1() );	
	double lat2 = xToLat( edge->getX2() );	
	double lon2 = yToLon( edge->getY2() );	

	double deltaLat = lat2-lat1;
	double deltaLon = lon2-lon1;
	
	double deltaT = .1;
	double t = 0;
	// First bracket zeros.
	double leftLat = lat1;
	double leftLon = lon1;
	double leftDer = haversineDeriv(leftLat, leftLon, 
									deltaLat, deltaLon, 
									weatherLat, weatherLon);
	double rightLat;
	double rightLon;
	double rightDer;
	for( t = 0; t < 1; t += deltaT)
	{
		rightLat = leftLat + deltaT*deltaLat;
		rightLon = leftLon + deltaT*deltaLon;
		rightDer = haversineDeriv( rightLat, rightLon, 
								   deltaLat, deltaLon, 
								   weatherLat, weatherLon);
		if(leftDer*rightDer <= 0)
		{
			break;
		}
	}
	// If we have exceeded our search space, then one of the endpoints is closer;
	if(leftDer*rightDer > 0)
	{
		double headDist = haversineDistance(lat1,lon1,weatherLat,weatherLon);
		double tailDist = haversineDistance(lat2,lon2,weatherLat,weatherLon);
		return std::min(headDist, tailDist);
	}

	// With a bracketed t in hand, we do a binary search
	double centerLat = (leftLat+rightLat)/2;
	double centerLon = (leftLon+rightLon)/2;


	while(abs(leftDer - rightDer) > EPSILON )
	{
		centerLat = (leftLat+rightLat)/2;
		centerLon = (leftLon+rightLon)/2;
		double centerDer = haversineDeriv(	centerLat, centerLon, 
											deltaLat, deltaLon, 
						  					weatherLat, weatherLon);
		// If zero is between center and left, move right to center
		if(centerDer*leftDer < 0) 
		{
			rightLat = centerLat;
			rightLon = centerLon;
			rightDer = centerDer;
		}
		else // move left point to center
		{
			leftLat = centerLat;
			leftLon = centerLon;
			leftDer = centerDer;
		}
	}

	return(haversineDistance(centerLat,centerLon,weatherLat,weatherLon) );

}

/**
	\brief Upper bound on the haversine distance from a point to a line

	We use this to save ourselves the time of numerically fiding the 
	point of closest approach between the weather and the edge.

	If the edge is (AC) and the weather point is (B) then
	there exists a t in [0,1] that minimizes d(B, A+t(A-C) )

	We denote by T the point A + t(A-C)

	Then d(B,T) \leq min( d(A,B) + d(A,T), d(C,B) + d(C,T) )

	Because t<.5, iff B is closer to A than C
	the minimum of d(A,T) and d(C,T) are less than 1/2 of d(A,C)

	d(B,T) \leq min( d(A,B) + d(A,C)/2, d(C,B) + d(A,C)/2 )

	\param edge The edge we're gauging distance from
	\param weatherLat Latitude of current weather
	\param weatherLon Longitude of current weather
*/
double TreeVerifier::boundDistWeatherToEdge( Edge *edge, double weatherLat, double weatherLon ) const
{
	double latA = xToLat( edge->getX1() );	
	double lonA = yToLon( edge->getY1() );	
	double latC = xToLat( edge->getX2() );	
	double lonC = yToLon( edge->getY2() );	

	double dAC = haversineDistance(latA, lonA, latC, lonC);
	double dAWeather  = haversineDistance(latA, lonA, weatherLat, weatherLon);
	double dCWeather  = haversineDistance(latC, lonC, weatherLat, weatherLon);

	return std::min( dAWeather, dCWeather) + dAC/2;
	
}

/**
	\brief Converts screen pixel x to latitude

	\param x Screen pixel's x coordinate
*/
double TreeVerifier::xToLat( double x ) const
{
	return centerLat + x*latPerPixel;
}

/**
	\brief Converts screen pixel y to longitude

	\param y Screen pixel's y coordinate
*/
double TreeVerifier::yToLon( double y ) const
{
	return centerLon + y*lonPerPixel;
}
