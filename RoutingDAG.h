// the class RoutingDAG is the graph on which we do the routing

#ifndef ROUTINGDAG_H
#define ROUTINGDAG_H

#define PI 3.14159265

#define TREE_NOT_GENERATED 1
#define TREE_GENERATED 2

#define NODES_READ_IN 1
#define NODES_NOT_READ_IN 2

#define EDGES_GENERATED 1
#define EDGES_NOT_GENERATED 2

#define BOTTOM_TREE 1
#define TAUTENED_TREE 2

// precomputed scale factor, length 1 in opengl coordinate system corresponds to 1.81081 nautical miles
#define NMILESPERPIXEL 1.81081

#include <vector>
#include "NodeAndEdge.h"

using namespace std;

class RoutingDAG
{
public:
	RoutingDAG();
	~RoutingDAG();
	// the key function, routing algorithm, used to generate a tree from a prebuilt graph
	bool generateTree(const vector<WeatherData*> &wDatas, vector<float> rnp, float effectiveThres, float routingThres);
	bool generateTautenedTree(const vector<WeatherData*> &wDatas, vector<float> rnp, float effectiveThres, float routingThres);
	void resetTree();
	void reset();
	void setminimumDistanceBetweenMergingNodes(double dis);
	void insertNode(Node* temp);
	void insertEdge(Edge* temp);
	void setNumLayers(int n);
	void setNodesReadInStatus(int status);
	int  getStatus();
	bool generateEdgeSet();
	void generateOperFlexPairs(float *radii, int length, vector<WeatherData*> &wDatas, float effectiveThres);
public:
	bool outputTreeInformation(double centerLati, double centerLong, double latiPerPixel, double longPerPixel, string &startTime, string &endTime);
private:
	int numLayers;
	// the entry nodes have to be in the order aligned along the outer boundary of the quadrant, with enough RNP separating them
	vector<Node*> entries;
	vector<Node*> nodes;
	vector<Node*> fixes;
	vector<Edge*> edges;
	// one number for each layer in the DAG, denoting which nodes are unavailable (Bottom most Filling), initialized to -1
	vector<int> layerStartingIndex;					// in the nodes vector, this helper vector records the starting position of each new layer 
	vector<int> layerUsedIndex;						// size is the number of layers
	vector<int> layerUsedIndexReverseDirection;		// when tautening the tree, we need one more usedIndex vector on the other side of the branch
private:
	int status;										// if a tree is generated
	int nodesReadIn;								// if the nodes of the DAG are ready
	int edgesGenerated;								// if we have generated all the edges
	int treeShapeStatus;							// if its bottommost tree or tautened tree
	double minimumDistanceBetweenMergingNodes;		// the minimum distance in between 2 merging nodes along a single branch
private:
	bool routeBranch(Node *start, int entryIndex, const vector<WeatherData*> &wDatas, float rnp, float effectiveThres, float routingThres);
	bool testRemainingBranchWhileMerging(Node *start, const vector<WeatherData*> &wDatas, float rnp, float effectiveThres, float routingThres);
	void setTreeBranchUp(Node* current, Node* start, float rnp);
	void setTreeBranchUpMerging(Node* start, float rnp);
	void treeBranchPostProcessing(Node* start, float rnp);
	bool testBranchWithEdge(Edge* current, float rnp, Node* start, int testType);
	bool testPreviousBranchEdge(Edge* current, float rnp, int entryIndex);
	bool testBranchWithNode(Node* current, float rnp, Node* start);
	bool testPreviousBranchNode(Node* current, float rnp, int entryIndex);
	bool testPreviousBranchTillCurrentLayerEdge(Edge* current, float rnp, int entryIndex);
	bool onPreviousBranch(Node* current, int entryIndex);
	bool onBranchStartingAt(Node* current, Node* start);
	bool generateLayerStartingIndexVector();
	int  findFeasiblePreviousEntryNode(int entryIndex);
	
	// helper functions for locating a node by its index or layer/layerIndex
	Node* fetchNode(int n);
	int getFetchNodeIndex(Node* temp);
	Node* findNode(int layer, int layerIndex);
	
	// helper functions when tautening the tree branches
	void updateLayerUsedIndexVector(int entryIndex);
	int findFeasibleNextEntryNode(int entryIndex);
	bool routeTautenedTreeBranch(Node *start, int entryIndex, const vector<WeatherData*> &wDatas, float rnp, 
								 float effectiveThres, float routingThres, int topMostTendency);
	bool onNextBranch(Node* current, int entryIndex);
	bool testNextBranchNode(Node* current, float rnp, int entryIndex);
	bool testNextBranchEdge(Edge* current, float rnp, int entryIndex);
	bool testnextBranchTillCurrentLayerEdge(Edge* current, float rnp, int entryIndex);
	void treeBranchPostProcessingForTreeTautening(Node* start, float rnp);
	void setTreeBranchUpForTreeTautening(Node* current, Node* start, float rnp);
	void setTreeBranchUpMergingForTreeTautening(Node *start, float rnp);
	Node* findMergingNodeOfPrevious2BranchesInBottomTree(int entryIndex);
	
	// helper functions used to test if the new branch satisfies the minimum distance between merging points constraint
	bool testDistanceTooCloseToMergingNodesOnNextBranch(Node* current, int entryIndex);
	bool testDistanceTooCloseToMergingNodesOnBranch(Node* current, Node* start);
	bool testDistanceTooCloseToMergingNodesOnPreviousBranch(Node* current, int entryIndex);
	bool testDistanceTooCloseToMergingNodesOnCurrentBranch(Node* toBeDecided, Node* current, int entryIndex);
	Node* nextNodeOfGivenNodeOnNextBranch(Node* current, int entryIndex);
	Node* nextNodeOfGivenNodeOnCurrentBranch(Node* current, int entryIndex);
};

#endif