#ifndef TREEVERIFIER_H
#define TREEVERIFIER_H

#include <string>
#include <vector>

#include "NodeAndEdge.h"
#include "WeatherData.h"
#include "RoutingDAG.h"
#include "UserInterface.h"

class TreeVerifier
{
private:
	const UserInterface &ui;
	double centerLat;
	double centerLon;
	double latPerPixel;
	double lonPerPixel;
public:
	TreeVerifier(const UserInterface &UI);

public: // verifier functions
	void appendReportToFiles(std::string invalidName, std::string validName);
	unsigned int countInvalidNodes() const;
	unsigned int countInvalidEdges() const;
private:
	bool doesNodeAvoidDangerousWeatherInEnsemble( Node *node, unsigned int ensembleIndex) const;
	bool doesEdgeAvoidDangerousWeatherInEnsemble( Edge *Edge, unsigned int ensembleIndex) const;

	double totalProbNodeIsSafe( Node *node ) const;
	double totalProbEdgeIsSafe( Edge *edge ) const;

	double shortestDistWeatherToEdge( Edge *edge, double weatherLat, double weatherLon ) const;
	double boundDistWeatherToEdge( Edge *edge, double weatherLat, double weatherLon ) const;

private: // Utility functions
	double xToLat( double x) const;
	double yToLon( double y) const;

};


#endif // already included

