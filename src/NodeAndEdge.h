/* define the Node class and Edge class, which are basic data structures that are used in the search graph*/

#ifndef NODEANDEDGE_H
#define NODEANDEDGE_H

#include <vector>
#include "WeatherData.h"

using namespace std;

#define ENTRY_NODE 1
#define INTERNAL_NODE 2
#define FIX_NODE 3

#define NOT_TREE_NODE 1
#define TREE_NODE 2

#define NOT_TREE_EDGE 1
#define TREE_EDGE 2

#define VISITED 1
#define NOT_VISITED 2

#define WEATHER_COLLISION 1
#define WEATHER_FREE 2
#define WEATHER_NOT_TESTED 3

class Edge;

/********************************************************************************************************************************************/
// define structure operationalflexity, which defines the probability for each radius of a center point or an edge

struct OperationalFlexibility
{
	double radius;
	double probability;
};

/********************************************************************************************************************************************/
// define class Node

class Node
{
public:
	Node(double ix, double iy, double iz, int itype = INTERNAL_NODE);
	~Node();
public:
	bool testRadiusWithWeatherDataSet(double r, const vector<WeatherData*> &wDatas, float effectiveThres, float routingThres);
	double testRadiusWithWeatherDataSet(double r, const vector<WeatherData*> &wDatas, float effectiveThres);
	bool collisionWithEdge(Edge *temp, float w);
	bool collisionWithNode(Node *temp, float w);
	void setTreeNode();				// set the node to be a node in the tree
	void setNotTreeNode();
	bool setLayer(int ilayer);
	bool setLayerIndex(int ilayerIndex);
	int getLayer();
	int getLayerIndex();
	bool treeNodeOrNot();
	
	// basic get and set functions
public:
	void insertInNodeEdge(Node* node, Edge* edge);						
	void insertOutNodeEdge(Node* node, Edge* edge);
	int getInSize();									// number of coming in edges
	int getOutSize();									// number of going out edges
	Node* getInNode(int index);
	Node* getOutNode(int index);
	Edge* getInEdge(int index);
	Edge* getOutEdge(int index);
	double getX();
	double getY();
	double getZ();
	void reset(double ix, double iy, double iz, int itype);
	void clearFreeRadiusVector();
	void setVisited(int n);
	int getVisited();
	bool addInDegree();
	int getInDegree();
	int getNodeType();
	void setDrawingRNP(float rnp);
	float getDrawingRNP();
	int getTreeOutEdgeIndex();
	void setTreeOutEdgeIndex(int index);
	void storeTreeOutEdgeIndexInformation();
	void restoreTreeOutEdgeIndexInformation();
	void resetTreeStatus();
	Node* getPrevTreeNode();
	Edge* getPrevTreeEdge();
	void setPrevTreeNode(Node* temp);
	void setPrevTreeEdge(Edge* temp);
	double getDistance();
	void setDistance(double dis);
	void resetWeatherCollisionStatus();
	void insertWeatherCollisionStatus(float rnp, int collisionStatus);
	int getWeatherCollisionStatus(float rnp);
	
	// Operational Flexibility Related Functions
public:
	int getFreeRadiusVecSize();
	bool getFreeRadiusResults(int n, double *radius, double *prob);	// get the value of nth element in the radius testing result (operational flexity)
	void insertFreeRadiusVec(float r, float prob);					// insert into the operational flexibility vector a pair of radius and probability
private:
	double x, y, z;										// the coordinates of the point
	vector<Node*> inNodes;								// the nodes that go into this node, size is the same as inEdges
	vector<Edge*> inEdges;								// the edges coming into this node, matching each of inNodes
	vector<Node*> outNodes;
	vector<Edge*> outEdges;								// the same as inNodes/inEdges, just the nodes/edges going out the this node
	int type;											// type 1, 2 or 3, entry, internal or fix
	int layer;											// which layer it is in
	int layerIndex;										// the index of the node inside its layer (from right to left in the quadrant)
	int treeOutEdgeIndex;								// which out edge is the edge in the tree if it is a tree node, which is also the index of out node
	int originalTreeOutEdgeIndex;						// used when tautening the tree, record the original bottommost tree information
	// data members that are used in the tree routing algorithm
private:
	int visited;										// mark if the node was visited already in the graph search algorithm
	int inDeg;											// currently, only allow in degree 2 trees
	int treeNode;										// after routing, if it is a node used in the tree
	float drawingRNP;									// the rnp when drawing, rnp is the width of the thick cylinder/2
	vector<OperationalFlexibility*> freeRadius;			// for a set of width, the probability that the node is free of obstacles
	Node* prevTreeNode;													
	Edge* prevTreeEdge;	
	double distance;									// when generating the tautened tree, used as distance in Dijkstra algorithm
	vector<float> weatherCollisionRNPs;					// define if a node is free of weather/ or collides with weather, record the tested rnps
	vector<int> weatherCollisionStatus;					// define if a node is free of weather/ or collides with weather for each float rnp
private:
	bool testRadiusWithWeatherData(double r, WeatherData* wData, float thres);	// test the surrounding of the point, if r=radius circle is weather free
	// test intersection between a disk whose center is (xC, yC), radius r, with a square, whose bottomleft corner is (x, y), side length c
	bool collisionBetweenDiskAndSquare(double xC, double yC, double r, double x, double y, double c);
	bool collisionWithEdgeHelper(float r, double ix1, double iy1, double iz1, double ix2, double iy2, double iz2, double iw);
	bool collisionWithNodeHelper(float r, double ix, double iy, double iz, double ir);
};

/********************************************************************************************************************************************/
// define class Edge

class Edge
{
public:
	Edge(Node* head, Node* tail);						// a directed edge, from node start to node end
	~Edge();
public:
	// algorithms related functions
	// test if an edge can be thicken enough to rnp*2 width and avoid weather obstacles
	bool testRNPWithWeatherDataSet(float rnp, const vector<WeatherData*> &wDatas, float effectiveThres, float routingThres);
	// variable "thres" in the following fuctions are used to denote which weathe data is considered to be hazardous
	// test the right side of the edge, if w=width rectangle is weather free
	bool testPathStretchWithWeatherDataSet(double width, const vector<WeatherData*> &wDatas, float effectiveThres, float routingThres);	
	// test the left side of the edge, if w=width rectangle is weather free
	bool testWiggleRoomWithWeatherDataSet(double width, const vector<WeatherData*> &wDatas, float effectiveThres, float routingThres);	
	
	// overloaded versions of the functions above, return the probability that an edge is clear
	double testRNPWithWeatherDataSet(float rnp, const vector<WeatherData*> &wDatas, float effectiveThres);
	double testPathStretchWithWeatherDataSet(double width, const vector<WeatherData*> &wDatas, float effectiveThres);
	double testWiggleRoomWithWeatherDataSet(double width, const vector<WeatherData*> &wDatas, float effectiveThres);
	
	bool collisionWithEdge(Edge *temp, float w);
	bool collisionWithNode(Node *temp, float w);
	void setTreeEdge();													// set the edge as an edge used in the tree
	void setNotTreeEdge();
	float getDrawingRNP();												
	void setDrawingRNP(float rnp);
public:
	// these functions are currently only used when outputting the path stretching and wiggle room information to an .xml file
	// compute the vertex coordinates of a rectangle on the right side of current edge, in clockwide order(head, tail, 3rd, 4th) 
	void computeRightSideRectangleVertices(double r, double* thirdX, double* thirdY, double* fourthX, double* fourthY);
	// compute the vertex coordinates of a rectangle on the left side of current edge, in clockwide order(1st, 2nd, tail, head) 
	void computeLeftSideRectangleVertices(double r, double* firstX, double* firstY, double* secondX, double* secondY);
public:
	// basic set and get functions
	Node* getHead();
	Node* getTail();
	double getLength();
	void reset();
	bool treeEdgeOrNot();
	void resetTreeStatus();
public:
	// information of RNP testing, path stretching testing and wiggle room testing results, oprational flexibility
	void insertOperFlex(float r, float prob, int vecIndex);
	void insertOperFlexDeviationCandidateNode(Node* temp);
	int getRNPVecSize();
	int getPathStretchingVecSize();
	int getWiggleRoomVecSize();
	int getDeviationNodesSize();
	bool getRNPResults(int n, double* r, double* prob);					// get the value of nth element in the rnpValues vector
	bool getPathStretchingResults(int n, double *width, double *prob);	// get the value of nth element in the pathStretching vector
	bool getWiggleRoomResults(int n, double *width, double *prob);	// get the value of nth element in the wiggleRoom vector
	Node* getDeviationNode(int n);
	void resetWeatherCollisionStatus();
	void insertWeatherCollisionStatus(float rnp, int collisionStatus);
	int getWeatherCollisionStatus(float rnp);
	/* set and get functions related to pathStretching and wiggleRoom defination here*/
private:
	double x1, y1, z1;									// the coordinates of its 2 endpoints
	double x2, y2, z2;
	Node* headNode;
	Node* tailNode;
	int layer, layerIndex;
	int treeEdge;
	float drawingRNP;									// the rnp that we are drawing, width/2
	double edgeLength;
	vector<OperationalFlexibility*> rnpValues;
	vector<OperationalFlexibility*> pathStretching;		// on the edge's right
	vector<OperationalFlexibility*> wiggleRoom;			// on the edge's left
	vector<Node*> deviationNodes;						// the nodes in this edge that are used to escape from the current tree
	// define if an edge is free of weather/ or collides with weather/ or not tested yet, avoid testing again(time consuming)
	vector<float> weatherCollisionRNPs;					// define if a node is free of weather/ or collides with weather, record the tested rnps
	vector<int> weatherCollisionStatus;					// define if a node is free of weather/ or collides with weather for each float rnp						
	
private:
	bool collisionWithWeatherDataHelper(float w, const vector<WeatherData*> &wDatas, float effectiveThres, float routingThres, int testType);
	double collisionWithWeatherDataHelper(float w, const vector<WeatherData*> &wDatas, float effectiveThres, int testType);
	
	// test the intersection between rectangle: center segments(x1, y1)->(x2, y2), width w*2, with square: bottomleft corner (x, y), side length c
	bool collisionBetweenRectangleAndSquare(double x1, double y1, double x2, double y2, double w, double x, double y, double c);
	bool collisionTestingHelper(double width, WeatherData *wData, float thres, int testType);
	bool collisionWithEdgeHelper(float w, double ix1, double iy1, double iz1, double ix2, double iy2, double iz2, double iw);
	bool collisionWithNodeHelper(float w, double ix, double iy, double iz, double ir);
};

#endif