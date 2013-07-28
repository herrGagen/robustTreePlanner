// the class RoutingDAG is the graph on which we do the routing
#ifndef ROUTINGDAG_H
#define ROUTINGDAG_H

#include <string>

#define PI 3.14159265

#define TREE_NOT_GENERATED 1
#define TREE_GENERATED 2

#define NODES_READ_IN 1
#define NODES_NOT_READ_IN 2

#define EDGES_GENERATED 1
#define EDGES_NOT_GENERATED 2

// #define OUTPUT_WEATHER_IN_XML

#define BOTTOM_TREE 1
#define TAUTENED_TREE 2

// precomputed scale factor, length 1 in opengl coordinate system corresponds to 1.81081 nautical miles
#define NMILESPERPIXEL 1.81081

#include <vector>
#include "NodeAndEdge.h"

class RoutingDAG
{
public:
	RoutingDAG();
	~RoutingDAG();
	// the key function, routing algorithm, used to generate 
        // a tree from a prebuilt graph
	bool generateTree(const std::vector<WeatherData> &wDataSets, std::vector<double> rnp, double effectiveThres, double routingThres);
	bool generateTautenedTree(const std::vector<WeatherData> &wDataSets, std::vector<double> rnp, double effectiveThres, double routingThres);
	void resetTree();
	void reset();
	void setminimumDistanceBetweenMergingNodes(double dis);
	void insertNode(Node* temp);
	void insertEdge(Edge* temp);
	void setNumLayers(unsigned int n);
	void setNodesReadInStatus(int status);
	int  getStatus();
	bool generateEdgeSet();
	void generateOperFlexPairs(const std::vector<double> &radii, const std::vector<WeatherData> &wDataSets, double effectiveThres);
public:
	bool outputTreeInformation(double centerLati, double centerLong, double latiPerPixel, double longPerPixel, const std::string &startTime, const std::string &endTime, const std::string &outputName, const std::vector<WeatherData> &wDataSets, double routingThresh);
private:
	unsigned int numLayers;
	// the entry nodes have to be in the order aligned along the outer boundary of the quadrant, with enough RNP separating them
	std::vector<Node*> entries;
	std::vector<Node*> nodes;
	std::vector<Node*> fixes;
	std::vector<Edge*> edges;
	// one number for each layer in the DAG, denoting which nodes are unavailable (Bottom most Filling), initialized to -1
	std::vector<int> layerStartingIndex;					// in the nodes std::vector, this helper std::vector records the starting position of each new layer 
	std::vector<int> layerUsedIndex;						// size is the number of layers
	std::vector<int> layerUsedIndexReverseDirection;		// when tautening the tree, we need one more usedIndex std::vector on the other side of the branch
private:
	int status;										// if a tree is generated
	int nodesReadIn;								// if the nodes of the DAG are ready
	int edgesGenerated;								// if we have generated all the edges
	int treeShapeStatus;							// if its bottommost tree or tautened tree
	double minimumDistanceBetweenMergingNodes;		// the minimum distance in between 2 merging nodes along a single branch
private:
	bool routeBranch(Node *start, unsigned int entryIndex, const std::vector<WeatherData> &wDataSets, double rnp, double effectiveThres, double routingThres);
	bool testRemainingBranchWhileMerging(Node *start, const std::vector<WeatherData> &wDataSets, double rnp, double effectiveThres, double routingThres);
	void setTreeBranchUp(Node* current, Node* start, double rnp);
	void setTreeBranchUpMerging(Node* start, double rnp);
	void treeBranchPostProcessing(Node* start, double rnp);
	bool testBranchWithEdge(Edge* current, double rnp, Node* start, int testType);
	bool testPreviousBranchEdge(Edge* current, double rnp, unsigned int entryIndex);
	bool testBranchWithNode(Node* current, double rnp, Node* start);
	bool testPreviousBranchNode(Node* current, double rnp, unsigned int entryIndex);
	bool testPreviousBranchTillCurrentLayerEdge(Edge* current, double rnp, unsigned int entryIndex);
	bool onPreviousBranch(Node* current, unsigned int entryIndex);
	bool onBranchStartingAt(Node* current, Node* start);
	bool generateLayerStartingIndexVector();
	int  findFeasiblePreviousEntryNode(unsigned int entryIndex);
	
	// functions for locating a node by its index or layer/layerIndex
	int getNodePointerIndex(Node* temp) const;
	Node* findNode(int layer, int layerIndex) const;

public:
	unsigned int getNumNodes() const;
	Node* getNodePointer(int n) const;
	unsigned int getNumEdges() const;
	Edge* getEdgePointer(int n) const;

private:
	// functions when tautening the tree branches
	void updateLayerUsedIndexVector(unsigned int entryIndex);
	int findFeasibleNextEntryNode(unsigned int entryIndex);
	bool routeTautenedTreeBranch(Node *start, unsigned int entryIndex, const std::vector<WeatherData> &wDataSets, double rnp, 
								 double effectiveThres, double routingThres, int topMostTendency);
	bool onNextBranch(Node* current, unsigned int entryIndex);
	bool testNextBranchNode(Node* current, double rnp, unsigned int entryIndex);
	bool testNextBranchEdge(Edge* current, double rnp, unsigned int entryIndex);
	bool testnextBranchTillCurrentLayerEdge(Edge* current, double rnp, unsigned int entryIndex);
	void treeBranchPostProcessingForTreeTautening(Node* start, double rnp);
	void setTreeBranchUpForTreeTautening(Node* current, Node* start, double rnp);
	void setTreeBranchUpMergingForTreeTautening(Node *start, double rnp);
	Node* findMergingNodeOfPrevious2BranchesInBottomTree(unsigned int entryIndex);
	
	// helper functions used to test if the new branch satisfies the minimum distance between merging points constraint
	bool testDistanceTooCloseToMergingNodesOnNextBranch(Node* current, unsigned int entryIndex);
	bool testDistanceTooCloseToMergingNodesOnBranch(Node* current, Node* start);
	bool testDistanceTooCloseToMergingNodesOnPreviousBranch(Node* current, unsigned int entryIndex);
	bool testDistanceTooCloseToMergingNodesOnCurrentBranch(Node* toBeDecided, Node* current, unsigned int entryIndex);
	Node* nextNodeOfGivenNodeOnNextBranch(Node* current, unsigned int entryIndex);
	Node* nextNodeOfGivenNodeOnCurrentBranch(Node* current, unsigned int entryIndex);

public:
	bool areAllNodesFarFromWeather( const std::vector<WeatherData> &wDataSets, 
											double rad, 
											double effectiveThresh, 
											double routingThresh );
};

#endif
