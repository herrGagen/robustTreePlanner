#include "RoutingDAG.h"
#include "NodeAndEdge.h"
#include <cmath>
#include <iostream>
#include <algorithm>
/********************************************************************************************************************************************/
// class Node

Node::Node(double ix, double iy, double iz, int itype)
{
	reset(ix, iy, iz, itype);
}

// reset THIS node
void Node::reset(double ix, double iy, double iz, int itype)
{
	x = ix;
	y = iy;
	z = iz;
	if(itype<ENTRY_NODE || itype>FIX_NODE)
	{
		type = INTERNAL_NODE;
	}
	else 
	{
		type = itype;
	}
	inDeg = 0;
	// initialize the other data members
	treeNode = NOT_TREE_NODE;
	layer = layerIndex = -1;
	drawingRNP = -1;
	visited = NOT_VISITED;
	// which out edge is the edge in the tree if it is a tree node, which is also the index of out node
	treeOutEdgeIndex = -1;								
	prevTreeNode = NULL;
	prevTreeEdge = NULL;
	inNodes.clear();
	outNodes.clear();
	inEdges.clear();
	outEdges.clear();
	freeRadius.clear();
	resetWeatherCollisionStatus();
}

Node::~Node()
{
}

// clear the free radius(for operational flexity) std::vector of a node
void Node::clearFreeRadiusVector()
{
	freeRadius.clear();
}

/**
	\brief test if the radius r ball centered at current node conflict with the weather data set
	
	This actually only does 2d testing if their z coordinates collide with each other

	\param effectiveThres if a weather cell's deviation probability is below this value, it's going to be considered as NULL
	\param routingThres means that we compute the weighted total probability of the weathercells p1*(0 or 1) + p2*(0 or 1) +..., 

	\retval True if collision occurs
	*/
bool Node::isDangerousWeatherCloserThanRadiusR(double r, const std::vector<WeatherData> &wDataSets, double effectiveThres, double routingThres)
{
#if !defined(DISABLE_COLLISION_CACHING)  
	if(getWeatherCollisionStatus(r) == WEATHER_COLLISION)			// if was tested to be colliding with the weather
	{
		return true;
	}
	if(getWeatherCollisionStatus(r) == WEATHER_FREE)				// no collision with weather
	{
		return false;										
	}
#endif
	if(wDataSets.size()==0)									// there is no weather read in
	{
		return false;										// the ball is considered to be clear of weather obstacle
	}
	double finalProbability = probabilityThatNodeIsClear(r,wDataSets, effectiveThres);

	if(finalProbability < routingThres)						// the probability that the weather is clear is not large enough
	{
		insertWeatherCollisionStatus(r, WEATHER_COLLISION);			// there is collision with the current weather data
		return true;										// the weather is too severe
	}
	insertWeatherCollisionStatus(r, WEATHER_FREE);
	return false;											// the node is considered to be clear of weather obstacle
}

// return the probability that the cell is clear
double Node::probabilityThatNodeIsClear(double r, const std::vector<WeatherData> &wDataSets, double effectiveThres)
{
	if(wDataSets.size()==0)									// there is no weather read in
	{
		return 1;											// there is no weather, therefore, there is 1 prabability that it's clear
	}
	double finalProbability = 0;
	for(unsigned int i=0; i<wDataSets.size(); i++)						// test the weather data members one by one
	{

		// if there is no intersection, means this weather data is "clear"
		if(!collidesWithWeatherDataSet(r, wDataSets[i], effectiveThres))
		{
			finalProbability += wDataSets[i].getProbability();					// add 1*probability of the weather data
		}
	}
	return finalProbability;
}


/**
\brief test if the radius r ball centered at current node conflict with the weatherdata

  This only does 2d testing if their z coordinates collide with each other

  This routine has been updated from treating weather data as the lower left corners
  of cellWidth*cellWidth squares.  Now they are circles of radius cellWidth.

  This returns true if any weather with a deviationProbability >= thresh is closer
  than r to any points in the weather ensemble wData.

  \param r Radius of node i.e. how far node must be from edge of weather.
  \param WeatherData A single weather ensemble.
  \param thresh Minimum deviationProbability for a weather cell to be dangerous
			

*/
bool Node::collidesWithWeatherDataSet(double r, const WeatherData &wData, double thres)
{
	if(wData.size() == 0)
		{
			return false;					// if the weather data does not exist, then just return false
		}
	double x1;
	double y1;
	double z1;
	double cellWidth;
	double cellHeight;
	double deviationProbability;
	for(unsigned int i=0; i<wData.size() ; i++)
	{
		// if successfully read out all the cell data, the bottomleft corner of the cell is (x1, y1, z1)
		if(wData.getCellData(i, &x1, &y1, &z1, &deviationProbability, &cellWidth, &cellHeight))	
		{
			if(deviationProbability <= thres)		// if the cell is not severe enough, simply skip
			{
				continue;						// test the next one							
			}
			if(z < z1 || z > z1+cellHeight)		// if the ranges in z direction have no overlap, then skip
			{
				continue;
			}
			else
			{
				double distanceInPixels = sqrt( (x-x1)*(x-x1) + (y-y1)*(y-y1) );
				double distanceInNM = distanceInPixels*NMILESPERPIXEL;
				// If it intersects with the current cell
				if(distanceInNM < cellWidth+r )	
				{
					return true;	
				}
			}
		}
	}
	return false;													// all the cells are tested and free of intersection
}

// reset all the elements related to the tree of a Node
void Node::resetTreeStatus()
{
	inDeg = 0;
	treeNode = NOT_TREE_NODE;
	drawingRNP = -1;
	freeRadius.clear();
	visited = NOT_VISITED;
}

// set THIS node to be a node of the tree
void Node::setTreeNode()
{
	treeNode = TREE_NODE;
}

void Node::setNotTreeNode()
{
	treeNode = NOT_TREE_NODE;
}

double Node::getDistance()
{	
	return distance;
}

void Node::setDistance(double dis)
{
	if(dis>=0)
	{
		distance = dis;
	}
	else
	{
		distance = 0;
	}
}

// reset the collision status to be never tested
void Node::resetWeatherCollisionStatus()
{
	weatherCollisionRNPs.clear();
	weatherCollisionStatus.clear();
}

// after a new testing, save the test results
void Node::insertWeatherCollisionStatus(double rnp, int collisionStatus)
{
	weatherCollisionRNPs.push_back(rnp);
	weatherCollisionStatus.push_back(collisionStatus);
}

int Node::getWeatherCollisionStatus(double rnp)
{
	for(unsigned int i=0; i<weatherCollisionRNPs.size(); i++)
	{
		if(weatherCollisionRNPs[i]>=rnp && weatherCollisionStatus[i]==WEATHER_FREE)
		{
			return WEATHER_FREE;
		}
		if(weatherCollisionRNPs[i]<=rnp && weatherCollisionStatus[i]==WEATHER_COLLISION)
		{
			return WEATHER_COLLISION;
		}
		if(weatherCollisionRNPs[i]==rnp)
		{
			return weatherCollisionStatus[i];
		}
	}
	return WEATHER_NOT_TESTED;				// didn't find the status for the specific rnp, meaning never tested before
}

// set the layer of the node in the search DAG
bool Node::setLayer(int ilayer)
{
	if(ilayer<0)
	{
		return false;
	}
	layer = ilayer;
	return true;
}

// set the layerIndex of the node in its layer
bool Node::setLayerIndex(int ilayerIndex)
{
	if(ilayerIndex<0)
	{
		return false;
	}
	layerIndex = ilayerIndex;
	return true;
}

int Node::getLayer()
{
	return layer;
}

int Node::getLayerIndex()
{
	return layerIndex;
}

// test if the node is a node in the tree
bool Node::isTreeNode()
{
	if(treeNode == TREE_NODE)
	{
		return true;
	}
	return false;
}

// only used when building the tree, served as prev[] in depth first search
void Node::setPrevTreeEdge(Edge *temp)
{
	prevTreeEdge = temp;
}

void Node::setPrevTreeNode(Node* temp)
{
	prevTreeNode = temp;
}

Node* Node::getPrevTreeNode()
{
	return prevTreeNode;
}

Edge* Node::getPrevTreeEdge()
{
	return prevTreeEdge;
}

// insert a coming in node and a coming in edge 
void Node::insertInNodeEdge(Node *node, Edge *edge)
{
	if(node && edge)			// if they are not NULL
	{
		inNodes.push_back(node);
		inEdges.push_back(edge);
	}
}

// insert a going out node and a going out edge 
void Node::insertOutNodeEdge(Node *node, Edge *edge)
{
	if(node && edge)
	{
		outNodes.push_back(node);
		outEdges.push_back(edge);
	}
}

unsigned int Node::getInSize()
{
	return inNodes.size();
}

unsigned int Node::getOutSize()
{
	return outNodes.size();
}

// search for the index-th node in the nodes that come into the node
Node* Node::getInNode(unsigned int index)
{
	if(index<getInSize() && index>=0)
	{
		return inNodes[index];
	}
	return NULL;
}

Node* Node::getOutNode(unsigned int index)
{
	if(index<getOutSize() && index>=0)
	{
		return outNodes[index];
	}
	return NULL;
}

Edge* Node::getInEdge(unsigned int index)
{
	if(index<getInSize() && index>=0)
	{
		return inEdges[index];
	}
	return NULL;
}

Edge* Node::getOutEdge(unsigned int index)
{
	if(index<getOutSize() && index>=0)
	{
		return outEdges[index];
	}
	return NULL;
}

double Node::getX()
{
	return x;
}

double Node::getY()
{
	return y;
}

double Node::getZ()
{
	return z;
}

bool Node::addInDegree()
{
	if(inDeg>=2)				// maximum indegree is 2
	{
		return false;			// cannot be added any more
	}
	inDeg++;
	return true;
}

int Node::getInDegree()
{
	return inDeg;
}

int Node::getNodeType()
{
	return type;
}

double Node::getDrawingRNP()
{
	return drawingRNP;
}

void Node::setDrawingRNP(double rnp)
{
	if(rnp>=0)
	{
		drawingRNP = rnp;
	}
}

void Node::setVisited(int n)
{
	if(n>=1 && n<=2)
	{
		visited = n;
	}
}

int Node::getVisited()
{
	return visited;
}

int Node::getTreeOutEdgeIndex()
{
	return treeOutEdgeIndex;
}

//only used when save the bottommost tree information by storing the original tree information in originalTreeOutEdgeIndex variable
void Node::storeTreeOutEdgeIndexInformation()
{
	originalTreeOutEdgeIndex = treeOutEdgeIndex;
}

//only used when restore the bottommost tree information by restoring the original tree information in originalTreeOutEdgeIndex variable
void Node::restoreTreeOutEdgeIndexInformation()
{
	treeOutEdgeIndex = originalTreeOutEdgeIndex;
}

// functions related to the free radius test, test if a specific radius is free of weather obstacles
int Node::getFreeRadiusVecSize()
{
	return freeRadius.size();
}

// functions related to the free radius test, test if a specific radius is free of weather obstacles
bool Node::getFreeRadiusResults(int n, double* radius, double* prob)
{
	if(n<0 || n>=getFreeRadiusVecSize())			// not that many std::vector elements, out of range
	{
		return false;
	}
	*radius = freeRadius[n]->radius;
	*prob = freeRadius[n]->probability;
	return true;
}

// insert into the operational flexibility std::vector a pair of radius and probability
void Node::insertFreeRadiusVec(double r, double prob)
{
	OperationalFlexibility* temp = new OperationalFlexibility;
	temp->radius = r;
	temp->probability = prob;
	freeRadius.push_back(temp);
}

void Node::setTreeOutEdgeIndex(unsigned int index)
{
	if(index<getOutSize() && index>=0)
	{
		treeOutEdgeIndex = index;
	}
}

// test if THIS node collides with another edge
bool Node::collisionWithEdge(Edge *temp, double r)
{
	return collisionWithEdgeHelper(r, temp->getHead()->getX(), temp->getHead()->getY(), temp->getHead()->getZ(), 
		temp->getTail()->getX(), temp->getTail()->getY(), temp->getTail()->getZ(), temp->getDrawingRNP());  
}

// test if THIS node collides with another edge
bool Node::collisionWithNode(Node* temp, double r)
{
	return collisionWithNodeHelper(r, temp->getX(), temp->getY(), temp->getZ(), temp->getDrawingRNP());
}

// test if THIS node(radius r) collides with another node, (ix, iy, iz), radius ir
bool Node::collisionWithNodeHelper(double r, double ix, double iy, double iz, double ir)
{
	if((x-ix)*(x-ix)+(y-iy)*(y-iy)+(z-iz)*(z-iz)>=(r+ir)*(r+ir))	// distance between the centers is greater than the sum of the radii
	{
		return false;				// no collision commited
	}
	return true;
}

// test if THIS node(radius r) collides with an edge, (ix1, iy1, iz2)->(ix2, iy2, iz2), width iw
// it's the same as testing the point (x, y, z) with the 3d segment (ix1, iy1, iz2)->(ix2, iy2, iz2), width iw+r
bool Node::collisionWithEdgeHelper(double r, double ix1, double iy1, double iz1, double ix2, double iy2, double iz2, double iw)
{
	double tWidth = r+iw;			// define the total width
	// first use bounding box to eleminate obvious cases
	if(std::min(iz1, iz2)-tWidth>=z || 
		std::max(iz1, iz2)+tWidth<=z || 
		std::min(iy1, iy2)-tWidth>=y || 
		std::max(iy1, iy2)+tWidth<=y ||
		std::min(ix1, ix2)-tWidth>=x || 
		std::max(ix1, ix2)+tWidth<=x)
	{
		return false;				// bounding box does NOT intersect the point
	}
	// check the intersection of the line passing (x, y, z) and is perpendicular to (ix1, iy1, iz2)->(ix2, iy2, iz2), test if its value is between 0 and 1
	double distanceSquare = (ix2-ix1)*(ix2-ix1)+(iy2-iy1)*(iy2-iy1)+(iz2-iz1)*(iz2-iz1);
	// compute the distance from (x, y, z) to the line defined by (ix1, iy1, iz2) and (ix2, iy2, iz2) and compare it to tWidth
	// tempT is used as the distance d here
	double vec[3] = {(y-iy1)*(z-iz2)-(z-iz1)*(y-iy2), (z-iz1)*(x-ix2)-(x-ix1)*(z-iz2), (x-ix1)*(y-iy2)-(y-iy1)*(x-ix2)};
	if (sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])/sqrt(distanceSquare) >= tWidth)			// 3d distance from a point to a line
	{
		return false;																			// distance is larger
	}
	// test if intersection is outside the range, AND the point is not too close to either end of the edge
	double tempT = -((ix1-x)*(ix2-ix1)+(iy1-y)*(iy2-iy1)+(iz1-z)*(iz2-iz1))/distanceSquare;
	if((tempT>=1 || tempT<=0) && (sqrt((x-ix1)*(x-ix1)+(y-iy1)*(y-iy1)+(z-iz1)*(z-iz1))>=tWidth) && (sqrt((x-ix2)*(x-ix2)+(y-iy2)*(y-iy2)+(z-iz2)*(z-iz2))>=tWidth))
	{
		return false;				// no way to intersect each other if the intersection point lies outside the range of the segment
	}

	return true;
}

/********************************************************************************************************************************************/
// class Edge

Edge::Edge(Node* head, Node* tail)
{
	if(head && tail)			// the pointers are non-NULL
	{
		headNode = head;
		tailNode = tail;
		x1 = head->getX();
		y1 = head->getY();
		z1 = head->getZ();
		x2 = tail->getX();
		y2 = tail->getY();
		z2 = tail->getZ();
		edgeLength = sqrt((z2-z1)*(z2-z1)+(y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
		layer = tailNode->getLayer();
		layerIndex = tailNode->getLayerIndex();
	}
	resetWeatherCollisionStatus();
	treeEdge = NOT_TREE_EDGE;
	drawingRNP= -1;
}

Edge::~Edge()
{
}

Node* Edge::getHead()
{
	return headNode;
}

Node* Edge::getTail()
{
	return tailNode;
}	

bool Edge::isTreeEdge()
{
	if(treeEdge==TREE_EDGE)
	{
		return true;
	}
	return false;
}

// set THIS edge as a part of the tree
void Edge::setTreeEdge()
{
	treeEdge = TREE_EDGE;
}

void Edge::setNotTreeEdge()
{
	treeEdge = NOT_TREE_EDGE;
}

double Edge::getLength()
{
	return edgeLength;
}

void Edge::reset()
{
	rnpValues.clear();
	pathStretching.clear();
	wiggleRoom.clear();
	deviationNodes.clear();
}

void Edge::resetTreeStatus()
{
	treeEdge = NOT_TREE_EDGE;
	drawingRNP = -1;								// the rnp that we are drawing, width/2
	rnpValues.clear();
	pathStretching.clear();							// on the edge's right
	wiggleRoom.clear();								// on the edge's left
	deviationNodes.clear();
}

double Edge::getDrawingRNP()
{
	return drawingRNP;
}

void Edge::setDrawingRNP(double rnp)
{
	if(rnp>0)
	{
		drawingRNP = rnp;
	}
}

// reset the weather collision status to never tested before
void Edge::resetWeatherCollisionStatus()
{
	weatherCollisionRNPs.clear();
	weatherCollisionStatus.clear();
}

// after a new testing, save the test results
void Edge::insertWeatherCollisionStatus(double rnp, int collisionStatus)
{
	weatherCollisionRNPs.push_back(rnp);
	weatherCollisionStatus.push_back(collisionStatus);
}

int Edge::getWeatherCollisionStatus(double rnp)
{
	for(unsigned int i=0; i<weatherCollisionRNPs.size(); i++)
	{
		if(weatherCollisionRNPs[i]>=rnp && weatherCollisionStatus[i]==WEATHER_FREE)
		{
			return WEATHER_FREE;
		}
		if(weatherCollisionRNPs[i]<=rnp && weatherCollisionStatus[i]==WEATHER_COLLISION)
		{
			return WEATHER_COLLISION;
		}
		if(weatherCollisionRNPs[i]==rnp)
		{
			return weatherCollisionStatus[i];
		}
	}
	return WEATHER_NOT_TESTED;				// didn't find the status for the specific rnp, meaning never tested before
}

// for each edge, we are going to compute a set of deviation points which can lead aircraft to leave the original tree
// this function returns the total number of such nodes
int Edge::getDeviationNodesSize()
{
	return deviationNodes.size();
}

// for each edge, we are going to compute a set of deviation points which can lead aircraft to leave the original tree
// this function returns a specific node based on its index
Node* Edge::getDeviationNode(int n)
{
	if(n<0 || n>=getDeviationNodesSize())
	{
		return NULL;
	}
	return deviationNodes[n];
}

// for each edge, we are going to compute a set of operational flexity pairs, for each width, the probability that the edge is free of weather obstacle
// this function returns the number of such pairs that we have, path stretching airspace is on the right of the current edge
int Edge::getPathStretchingVecSize()
{
	return pathStretching.size();
}

// for each edge, we are going to compute a set of operational flexity pairs, for each width, the probability that the edge is free of weather obstacle
// this function returns such a pair based on its index
bool Edge::getPathStretchingResults(int n, double *width, double *prob)
{
	if(n<0 || n>getPathStretchingVecSize())
	{
		return false;
	}
	*prob = pathStretching[n]->probability;
	*width = pathStretching[n]->radius;
	return true;
}

// for each edge, we are going to compute a set of operational flexity pairs, for each width, the probability that the edge is free of weather obstacle
// this function returns the number of such pairs that we have, wiggle room airspace is on the left of the current edge
int Edge::getWiggleRoomVecSize()
{
	return wiggleRoom.size();
}

// for each edge, we are going to compute a set of operational flexity pairs, for each width, the probability that the edge is free of weather obstacle
// this function returns such as pair based on its index
bool Edge::getWiggleRoomResults(int n, double* width, double* prob)
{
	if(n<0 || n>getWiggleRoomVecSize())
	{
		return false;
	}
	*prob = wiggleRoom[n]->probability;
	*width = wiggleRoom[n]->radius;
	return true;
}

// for each edge, we are going to compute a set of operational flexity pairs, for each width, the probability that the edge is free of weather obstacle
// this function returns the number of such pairs that we have, rnp values are on both sides of the edge
int Edge::getRNPVecSize()
{
	return rnpValues.size();
}

// for each edge, we are going to compute a set of operational flexity pairs, for each width, the probability that the edge is free of weather obstacle
// this function returns such as pair based on its index
bool Edge::getRNPResults(int n, double *r, double *prob)
{
	if(n<0 || n>getRNPVecSize())
	{
		return false;
	}
	*prob = rnpValues[n]->probability;
	*r = rnpValues[n]->radius;
	return true;
}

// insert an operational flexibility pair into corresponding std::vector (marked by vecIndex value)
void Edge::insertOperFlex(double r, double prob, int vecIndex)
{
	if(vecIndex<1 || vecIndex>3)
	{
		std::cerr << "\nCannot insert the operational plexibility pair into a valid std::vector..."<<std::endl;
		return;
	}
	OperationalFlexibility* temp = new OperationalFlexibility;
	temp->probability = prob;
	temp->radius = r;
	switch (vecIndex)
	{
	case 1:
		rnpValues.push_back(temp);
		break;
	case 2:
		pathStretching.push_back(temp);
		break;
	case 3:
		wiggleRoom.push_back(temp);
		break;
	default:
		break;
	}
}

void Edge::insertOperFlexDeviationCandidateNode(Node *temp)
{
	deviationNodes.push_back(temp);
}

// test the rnp on both side of THIS edge with the weatherDataSet (a set of boxes)
bool Edge::isAnyWeatherWithinLaneWidthW(double rnp, const std::vector<WeatherData> &wDataSets, double effectiveThres, double routingThres)
{
#if !defined(DISABLE_COLLISION_CACHING)	
  if(getWeatherCollisionStatus(rnp) == WEATHER_COLLISION)			// if was tested to be colliding with the weather
	{
		return true;
	}
	if(getWeatherCollisionStatus(rnp) == WEATHER_FREE)				// no collision with weather
	{
		return false;
	}
#endif
	if(wDataSets.size()==0)	
	{
		return false;
	}
	if(collisionWithWeatherCheck(rnp, wDataSets, effectiveThres, routingThres, 1))
	{
		insertWeatherCollisionStatus(rnp, WEATHER_COLLISION);
		return true;
	}
	insertWeatherCollisionStatus(rnp, WEATHER_FREE);
	return false;
}
// test the right side of THIS edge with the WeatherDataSet (a set of boxes), the width of the rectangle is passed into the function as parameter width
bool Edge::testPathStretchWithWeatherDataSet(double width, const std::vector<WeatherData> &wDataSets, double effectiveThres, double routingThres)
{
	if(wDataSets.size()==0)	
	{
		return false;
	}
	return collisionWithWeatherCheck(width, wDataSets, effectiveThres, routingThres, 2);
}
// test the left side of THIS edge with the WeatherDataSet (a set of boxes)
bool Edge::testWiggleRoomWithWeatherDataSet(double width, const std::vector<WeatherData> &wDataSets, double effectiveThres, double routingThres)
{
	if(wDataSets.size()==0)
	{
		return false;
	}
	return collisionWithWeatherCheck(width, wDataSets, effectiveThres, routingThres, 3);
}

// overloaded function to test the probability that THIS edge is clear of obstacles
double Edge::sumOfProbOfAllWeatherCloserThanLaneWidthW(double rnp, const std::vector<WeatherData> &wDataSets, double effectiveThres)
{
	if(wDataSets.size()==0)	
	{
		return 1;
	}
	return probabilityThatEdgeIsClear(rnp, wDataSets, effectiveThres, 1);
}

// overloaded function to test the probability that the right side of THIS edge is clear of obstacles
double Edge::testPathStretchWithWeatherDataSet(double width, const std::vector<WeatherData> &wDataSets, double effectiveThres)
{
	if(wDataSets.size()==0)	
	{
		return 1;
	}
	return probabilityThatEdgeIsClear(width, wDataSets, effectiveThres, 2);
}

// overloaded function to test the probability that the left side of THIS edge is clear of obstacles
double Edge::testWiggleRoomWithWeatherDataSet(double width, const std::vector<WeatherData> &wDataSets, double effectiveThres)
{
	if(wDataSets.size()==0)
	{
		return 1;
	}
	return probabilityThatEdgeIsClear(width, wDataSets, effectiveThres, 3);
}

// test rnp, pathStretching or Wiggle airspace are basically the same, only difference is the rectangle that we try to test
// parameter testType means the type we are testing, 1 for rnp, 2 for pathStretching, 3 for Wiggle Airspace. 
// test the values with a set of weather data files stored in a std::vector of weatherdata pointers
bool Edge::collisionWithWeatherCheck(double w, const std::vector<WeatherData> &wDataSets, double effectiveThres, double routingThres, int testType)
{
	if(wDataSets.size()==0)										// no weather data to test
	{
		return false;											// no intersection with weather data
  }
	double finalProbability = probabilityThatEdgeIsClear(w, wDataSets, effectiveThres, testType);
	if(finalProbability < routingThres)							// compare the probability of clear weather with the routing required probability
	{
		return true;											// the weather is too severe
	}
	return false;												// the weather is clear
}

// Returns the final probability value (used for output)
double Edge::probabilityThatEdgeIsClear(double distFromEdge, const std::vector<WeatherData> &wDataSets, double effectiveThres, int testType)
{
	if(wDataSets.size()==0)										// no weather data to test
	{
		return 1;											// no intersection with weather data
	}
	double finalProbability = 0;							// the total probability of all impacted weather
	for(unsigned int i=0; i<wDataSets.size(); i++)
	{
		if( !collidesWithWeatherDataSet(distFromEdge, wDataSets[i], effectiveThres, testType) )	 // if there is no collision
		{
			finalProbability += wDataSets[i].getProbability();
		}
	}												// the weather is clear
	return finalProbability;
}

/**
   test rnp, pathStretching or Wiggle airspace are basically the same, only difference is the rectangle that we try to test
 
 \param testType means the type we are testing, 1 for rnp, 2 for pathStretching, 3 for Wiggle Airspace
*/
bool Edge::collidesWithWeatherDataSet(double width, const WeatherData &wSet, double thres, int testType)
{
	unsigned int numCells = wSet.size();
	double x;
	double y;
	double z;
	double cellWidth;
	double cellHeight;
	double deviationProbability;
	for(unsigned int i=0; i<numCells; i++)
	{
		// if successfully read out all the cell data, the bottomleft corner of the cell is (x, y, z)
		if(wSet.getCellData(i, &x, &y, &z, &deviationProbability, &cellWidth, &cellHeight))	
		{
			if(deviationProbability<=thres)		// if the cell is not severe enough, simply skip
			{
				continue;						// test the next one	
			}
			// the ranges in z direction have no overlap, then skip
			if(std::min(z1, z2)>=z+cellHeight || std::max(z1, z2)<=z)
			{
				continue;						
			}
			else
			{
				// there is overlap in z direction, therefore, test only the part of the edge that overlaps with the box in z direction
				double tempX1; 
				double tempY1;
				double tempX2;
				double tempY2;				// define the 2 endpoints of the part of the edge
				if(z1==z2)
				{
					z1 = z2+0.0001;							// z1 and z2 should be different, otherwise may cause confusion
				}
				// first, the vertex that is near the smaller z of the 2 points, test is z drops in between the 2 endpoints or below the lower of the 2 points
				if(z>std::min(z1, z2) && z<std::max(z1, z2))					// if z is inbetween the segment
				{
					tempX2 = x1+(x2-x1)*std::abs(z-z1)/std::abs(z2-z1);
					tempY2 = y1+(y2-y1)*std::abs(z-z1)/std::abs(z2-z1);
				}
				else												// meaning z is even smaller than std::min(z1, z2)
				{
					tempX2 = (std::min(z1, z2)==z1)? x1 : x2;
					tempY2 = (std::min(z1, z2)==z1)? y1 : y2;
				}
				// first, the vertex that is near the larger z of the 2 points
				if(z+cellHeight>std::min(z1, z2) && z+cellHeight<std::max(z1, z2))					// if z+cellHeight is inbetween the segment
				{
					tempX1 = x1+(x2-x1)*std::abs(z+cellHeight-z1)/std::abs(z2-z1);
					tempY1 = y1+(y2-y1)*std::abs(z+cellHeight-z1)/std::abs(z2-z1);
				}
				else												// meaning z is even greater than std::max(z1, z2)
				{
					tempX1 = (std::min(z1, z2)==z1)? x2 : x1;
					tempY1 = (std::min(z1, z2)==z1)? y2 : y1;
				}
        if(testType == 1)
				{
					// test RNP only: the order of the 2 points is not important
					// now the segment (tempX1, tempY1)<->(tempX2, tempY2) is the part of the original edge that we are interested in detecting collision

          if( collisionBetweenRectangleAndCircle(tempX1, tempY1, tempX2, tempY2, width, x, y, cellWidth))	// if it intersects with the current cell
					{
						return true;																// then there is intersection
					}
				}
				else
				{
					if(std::min(z1, z2)==z1)					// the edge must be from head to tail in order to talk about RIGHT or LEFT size of the edge
					{									// to keep the direction, swap the 2 endpoints of the part of THIS edge that we are interested in
						double swap1 = tempX1;			double swap2 = tempY1;
						tempX1 = tempX2;				tempY1 = tempY2;
						tempX2 = swap1;					tempY2 = swap2;
					}
					double vectorLength = sqrt((tempX1-tempX2)*(tempX1-tempX2)+(tempY1-tempY2)*(tempY1-tempY2));
					width = width/2;
					// test path stretching on the right side of THIS edge (x1, y1, z1) to (x2, y2, z2)
					// now test the right side of this edge, where the new center line is defined by points 
					// (x1 + (y2-y1)*width/vectorLength, y1+(x1-x2)*width/vectorLength) and (x2 + (y2-y1)*width/vectorLength, y2+(x1-x2)*width/vectorLength)
					if(testType == 2)
					{
						if(collisionBetweenRectangleAndCircle(tempX1+(tempY2-tempY1)*width/vectorLength, 
							tempY1+(tempX1-tempX2)*width/vectorLength, 
							tempX2+(tempY2-tempY1)*width/vectorLength, 
							tempY2+(tempX1-tempX2)*width/vectorLength, 
							width, 
							x, 
							y, 
							cellWidth))
						{
							return true;					// there is collision
						}
					}
					// test path stretching on the left side of THIS edge (x1, y1, z1) to (x2, y2, z2)
					// now test the left side of this edge, where the new center line is defined by points 
					// (x1 + (y1-y2)*width/vectorLength, y1+(x2-x1)*width/vectorLength) and (x2 + (y1-y2)*width/vectorLength, y2+(x2-x1)*width/vectorLength)
					if(testType == 3)
					{
						if(collisionBetweenRectangleAndCircle(tempX1+(tempY1-tempY2)*width/vectorLength, 
							tempY1+(tempX2-tempX1)*width/vectorLength, 
							tempX2+(tempY1-tempY2)*width/vectorLength, 
							tempY2+(tempX2-tempX1)*width/vectorLength, 
							width, 
							x, 
							y, 
							cellWidth) )
						{
							return true;					// there is collision
						}
					}
				}
			}
		}			// if after testing all the weather cells and the function didn't return, then it means no collision, return false
	}	
	return false;
}


/**
	\brief Tests if edge with with laneWidth w collides with weather at (x,y) with radius c.

	Note: w and c are treated as radii, not as diameters.  So if they need to be halved, this
	function, or whatever is calling it, needs to be changed.

	We test if an edge and weather collide by:

	1) Parameterizing the edge as the line segment: (ex(t), ey(t)) =  (x1,y1) + t( (x2-x1), (y2-y1) )

	2) Finding the t that minimizes the distance between that line segment and (x,y).  ( via calculus )

	3) If closestT  \in [0,1] we check this distance, otherwise, we check min of dista at t=0 and t=1

	4) Return true if distance from step 3 is less than allowed distance of w + c

	\param x1 
	\retval True if weather collides with edge.


*/
bool Edge::collisionBetweenRectangleAndCircle(double x1Pixels, double y1Pixels, double x2Pixels, double y2Pixels, double w, double xPixels, double yPixels, double c)
{
	double x1 = x1Pixels*NMILESPERPIXEL;
	double y1 = y1Pixels*NMILESPERPIXEL;
	double x2 = x2Pixels*NMILESPERPIXEL;
	double y2 = y2Pixels*NMILESPERPIXEL;
	double x  =  xPixels*NMILESPERPIXEL;
	double y  =  yPixels*NMILESPERPIXEL;
	

	double deltaX = x2-x1;
	double deltaY = y2-y1;

	/*
	Here we are
	1) Parameterizing the edge as the line segment: (ex(t), ey(t)) =  (x1,y1) + t( (x2-x1), (y2-y1) )
	2) Finding the t that minimizes the distance between that line segment and (x,y).  ( via calculus )

	The distance^2 between (x1,y1) + t( (x2-x1), (y2-y1) ) and (x,y) is:
	(x-x1-t*deltaX)^2 + (y - y1 - t*deltaY)^2

	The partial derivative of this wrt t is:
	-2*deltaX*t*(x-x1-t*deltaX) -2*deltaY*t*(y-y1-t*deltaY)

	Setting this to zero and solving for t:
	0 = -2*deltaX*t*(x-x1-t*deltaX) - 2*deltaY*t(y-y1-t*deltaY)
	0 = deltaX*(x-x1-t*deltaX) +deltaY*(y-y1-t*deltaY)
	t * (deltaX^2 + deltaY^2) = (x-x1)*deltaX + (y-y1)*deltaY;
	t = (  (x-x1)*deltaX + (y-y1)*deltaY ) / (deltaX^2 + deltaY^2);
	*/
	double closestT = (  (x-x1)*deltaX + (y-y1)*deltaY ) / (deltaX*deltaX + deltaY*deltaY);

	// 3) If closestT  \in [0,1] we check this distance, otherwise, we check min of dista at t=0 and t=1
	double closestDistance;
	if(closestT >=0 && closestT <= 1)
	{
		double xt = x1 + closestT*deltaX;
		double yt = y1 + closestT*deltaY;
		double distanceT = sqrt( (x-xt)*(x-xt) + (y-yt)*(y-yt) );
		closestDistance = distanceT;
	}
	else // if t=0, we are at (x1,y1) and if t=1 we are at (x2,y2)
	{
		double distance1 = sqrt( (x-x1)*(x-x1) + (y-y1)*(y-y1) );
		double distance2 = sqrt( (x-x2)*(x-x2) + (y-y2)*(y-y2) );
		closestDistance = std::min(distance1, distance2);
	}
	// Finally return if the distance between the weather point and edge is less than the allowed distance, 
	return (closestDistance < w + c);
}

// test if THIS edge collides with another edge
bool Edge::collisionWithEdge(Edge *temp, double w)
{
	return collisionWithEdgeHelper(w, temp->getHead()->getX(), temp->getHead()->getY(), temp->getHead()->getZ(), 
		temp->getTail()->getX(), temp->getTail()->getY(), temp->getTail()->getZ(), temp->getDrawingRNP()); 
}

// test if THIS edge collides with another node
bool Edge::collisionWithNode(Node *temp, double w)
{
	return collisionWithNodeHelper(w, temp->getX(), temp->getY(), temp->getZ(), temp->getDrawingRNP());
}

// test if THIS edge(width w) collides with another node, (ix, iy, iz), radius ir
bool Edge::collisionWithNodeHelper(double w, double ix, double iy, double iz, double ir)
{
	double tWidth = w+ir;			// define the total width
	// first use bounding box to eleminate obvious cases
	if( std::min(z1, z2)-tWidth>=iz || 
		std::max(z1, z2)+tWidth<=iz || 
		std::min(y1, y2)-tWidth>=iy || 
		std::max(y1, y2)+tWidth<=iy ||
		std::min(x1, x2)-tWidth>=ix || 
		std::max(x1, x2)+tWidth<=ix)
	{
		return false;				// bounding box does NOT intersect the point
	}
	// check the intersection of the line passing (x, y, z) and is perpendicular to (ix1, iy1, iz2)->(ix2, iy2, iz2), test if its value is between 0 and 1
	double distanceSquare = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
	// compute the distance from (x, y, z) to the line defined by (ix1, iy1, iz2) and (ix2, iy2, iz2) and compare it to tWidth
	double vec[3] = {(iy-y1)*(iz-z2)-(iz-z1)*(iy-y2), (iz-z1)*(ix-x2)-(ix-x1)*(iz-z2), (ix-x1)*(iy-y2)-(iy-y1)*(ix-x2)};
	if (sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])/sqrt(distanceSquare) >= tWidth)			// 3d distance from a point to a line
	{
		return false;																			// distance is larger
	}
	// test if intersection is outside the range, AND the point is not too close to either end of the edge
	double tempT = -((x1-ix)*(x2-x1)+(y1-iy)*(y2-y1)+(z1-iz)*(z2-z1))/distanceSquare;
	if( (tempT>=1 || tempT<=0) && 
		(sqrt((x1-ix)*(x1-ix)+(y1-iy)*(y1-iy)+(z1-iz)*(z1-iz))>=tWidth) && 
		(sqrt((x2-ix)*(x2-ix)+(y2-iy)*(y2-iy)+(z2-iz)*(z2-iz))>=tWidth))	
	{
		return false;				// no way to intersect each other if the intersection point lies outside the range of the segment
	}
	return true;
}

// test if THIS edge(width w) collides with another edge, (ix1, iy1, iz2)->(ix2, iy2, iz2), width iw
// test the collision of a cylinder with radii (w+iw) with the THIS edge
bool Edge::collisionWithEdgeHelper(double w, double ix1, double iy1, double iz1, double ix2, double iy2, double iz2, double iw)
{
	double tWidth = w+iw;
	if( std::min(ix1, ix2)-tWidth > std::max(x1, x2) || 
		std::max(ix1, ix2)+tWidth < std::min(x1, x2) || 
		std::min(iy1, iy2)-tWidth > std::max(y1, y2) || 
		std::max(iy1, iy2)+tWidth < std::min(y1, y2) ||
		std::min(iz1, iz2)-tWidth > std::max(z1, z2) || 
		std::max(iz1, iz2)+tWidth < std::min(z1, z2) )
	{
		return false;							// bounding box does NOT intersect the segment
	}
	// test the distance between 2 lines, distance formula between 2 lines
	double vec[3] = {(y2-y1)*(iz2-iz1)-(z2-z1)*(iy2-iy1), (z2-z1)*(ix2-ix1)-(x2-x1)*(iz2-iz1), (x2-x1)*(iy2-iy1)-(y2-y1)*(ix2-ix1)};
	// the distance is std::abs((ix1-x1)*vec[0]+(iy1-y1)*vec[1]+(iz1-z1)*vec[2])/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]
	if (std::abs((ix1-x1)*vec[0]+(iy1-y1)*vec[1]+(iz1-z1)*vec[2])/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])>=tWidth)
	{
		return false;
	}
	// distance is too large for a collision to exist, then what if the distance is not that large??? test if the intersection happens in the range of the 2 segments
	// 2 lines: L1 = X1+(X2-X1)*s, L2 = X3 + (X4-X3)*s are points on the lines, we find the common perpendicular by computing 2 perpendicular planes of 2 lines respectively,
	// then the other point must be on the plane.
	double temp1 = (ix2-ix1)*(x2-x1)+(iy2-iy1)*(y2-y1)+(iz2-iz1)*(z2-z1);
	double temp2 = (x2-x1)*(ix1-x1)+(y2-y1)*(iy1-y1)+(z2-z1)*(iz1-z1);
	double distanceSquare = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
	double r = ((ix2-ix1)*(x1-ix1)+(iy2-iy1)*(y1-iy1)+(iz2-iz1)*(z1-iz1)+temp1*temp2/distanceSquare)/
		((ix2-ix1)*(ix2-ix1)+(iy2-iy1)*(iy2-iy1)+(iz2-iz1)*(iz2-iz1)-temp1*temp1/distanceSquare);
	if(r>=1 || r<=0)					// the parameter of the intersection between the common perpendicular and the input segment
	{
		return false;					// if the intersection is not on the segment, then no intersection
	}
	double s = temp2/distanceSquare+temp1*r/distanceSquare;
	if(s>=1 || s<=0)					// the parameter of the intersection between the common perpendicular and THIS segment
	{
		return false;					// if the intersection is not on the segment, then no intersection
	}
	return true;						// else, must be intersected
}

// these functions are currently only used when outputting the path stretching and wiggle room information to an .xml file
// compute the vertex coordinates of a rectangle on the right side of current edge, in clockwide order(head, tail, 3rd, 4th), rectangle width is r
void Edge::computeRightSideRectangleVertices(double r, double* thirdX, double* thirdY, double* fourthX, double* fourthY)
{
	*thirdX = x2 + (y2-y1)*r/edgeLength;
	*thirdY = y2 + (x1-x2)*r/edgeLength;
	*fourthX= x1 + (y2-y1)*r/edgeLength;
	*fourthY = y1 + (x1-x2)*r/edgeLength;
}

// compute the vertex coordinates of a rectangle on the left side of current edge, in clockwide order(1st, 2nd, tail, head), rectangle width is r
void Edge::computeLeftSideRectangleVertices(double r, double* firstX, double* firstY, double* secondX, double* secondY)
{
	*secondX = x2 - (y2-y1)*r/edgeLength;
	*secondY = y2 - (x1-x2)*r/edgeLength;
	*firstX = x1 - (y2-y1)*r/edgeLength;
	*firstY = y1 - (x1-x2)*r/edgeLength;
}