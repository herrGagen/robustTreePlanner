#include "RoutingDAG.h"
#include <stack>
#include <deque>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "math.h"
#include <string>

RoutingDAG::RoutingDAG()
{
	// vectors are empty at the beginning
	status = TREE_NOT_GENERATED;
	nodesReadIn = NODES_NOT_READ_IN;
	edgesGenerated = EDGES_NOT_GENERATED;
	treeShapeStatus = BOTTOM_TREE;
	minimumDistanceBetweenMergingNodes=0;		// on default, there is no requirement on this
}

// delete all the memory that is allocated to the nodes and edges
RoutingDAG::~RoutingDAG()
{
	for(int i=0; i<entries.size(); i++)
		delete entries[i];
	for(int i=0; i<nodes.size(); i++)
		delete nodes[i];
	for(int i=0; i<fixes.size(); i++)
		delete fixes[i];
	for(int i=0; i<edges.size(); i++)
		delete edges[i];
	entries.clear();
	nodes.clear();
	fixes.clear();
	edges.clear();
	layerUsedIndex.clear();
	layerStartingIndex.clear();
}

// insert an edge into the graph
void RoutingDAG::insertEdge(Edge *temp)
{
	edges.push_back(temp);
}

// set the number of layers in the routingDAG
void RoutingDAG::setNumLayers(int n)
{
	if(n>=2)					// at least 2 levels of fix nodes and entry nodes
	{
		numLayers = n;
		for(int i=0; i<numLayers; i++)			// initialize the layerUsedIndex vector which will be used when computing the tree
			layerUsedIndex.push_back(-1);
		for(int i=0; i<fixes.size(); i++)
			fixes[i]->setLayer(n-1);			// all the fix nodes now know that they are in layer (n-1)
	}
}

void RoutingDAG::setNodesReadInStatus(int status)
{
	if(status>=1 && status<=2)
		nodesReadIn = status;
	else nodesReadIn = NODES_NOT_READ_IN;
}

// insert a node into the graph
void RoutingDAG::insertNode(Node *temp)
{
	switch(temp->getNodeType())
	{
		case INTERNAL_NODE:
			nodes.push_back(temp);
			break;
		case ENTRY_NODE:
			entries.push_back(temp);
			break;
		case FIX_NODE:
			fixes.push_back(temp);
			break;
		default:
			break;
	}
}

// reset every node and edge to be NOT in the tree
// reset the weather status to be WEATHER_NOT_TESTED, a helper scheme to reduce the number of testings against weather data
void RoutingDAG::resetTree()
{
	for(int i=0; i<entries.size(); i++)
	{
		entries[i]->resetTreeStatus();
		entries[i]->resetWeatherCollisionStatus();
	}
	for(int i=0; i<nodes.size(); i++)
	{
		nodes[i]->resetTreeStatus();
		nodes[i]->resetWeatherCollisionStatus();
	}
	for(int i=0; i<fixes.size(); i++)
	{
		fixes[i]->resetTreeStatus();
		fixes[i]->resetTreeStatus();
	}
	for(int i=0; i<edges.size(); i++)
	{
		edges[i]->resetTreeStatus();
		edges[i]->resetWeatherCollisionStatus();
	}
	for(int i=0; i<layerUsedIndex.size(); i++)
		layerUsedIndex[i] = -1;
	status = TREE_NOT_GENERATED;
}

void RoutingDAG::reset()
{
	for(int i=0; i<entries.size(); i++)
		delete entries[i];
	for(int i=0; i<nodes.size(); i++)
		delete nodes[i];
	for(int i=0; i<fixes.size(); i++)
		delete fixes[i];
	for(int i=0; i<edges.size(); i++)
		delete edges[i];
	entries.clear();
	nodes.clear();
	fixes.clear();
	edges.clear();
	layerUsedIndex.clear();
	layerStartingIndex.clear();
	status = TREE_NOT_GENERATED;
	nodesReadIn = NODES_NOT_READ_IN;
	edgesGenerated = EDGES_NOT_GENERATED;
	treeShapeStatus = BOTTOM_TREE;
}

int RoutingDAG::getStatus()						// get the tree generated status
{
	return status;
}

// set the minimumDistanceBetweenMergingNodes parameter, which is a constraint when generating trees
// the constraint: two merging node on the same branch have to be some distance apart from each other
// the input was in nmile unit, so we convert the length it to screen opengl coordinate system first
void RoutingDAG::setminimumDistanceBetweenMergingNodes(double dis)
{
	if(dis >= 0)
	{
		minimumDistanceBetweenMergingNodes = dis;
	}
}

// output the tree information into a .tre ASCII file
// parametres are used to convert screen coordinates to lati/longs, and the time range of testing is also passed in as parameters
bool RoutingDAG::outputTreeInformation(double centerLati, double centerLong, double latiPerPixel, double longPerPixel, string &startTime, string &endTime, string &outputName)
{
	if(status==TREE_NOT_GENERATED)
	{
		std::cerr << "\nPlease generate the tree first before outputting the tree information!" << std::endl;
		return false;
	}
	string XMLFileName("out.xml");
        ofstream os;
        do
          {
            string tempFileName;
            tempFileName = outputName.c_str();
            std::cout << outputName.c_str() << std::endl;
            // std::cout<<"\nPlease provide a name for the output file:\n";
            // std::getline(cin,tempFileName);
            if( tempFileName.length() > 0 )
              {
				size_t xmlLocation;

				// Append extension if required
				xmlLocation = tempFileName.find(".xml");

				// Not found, try all caps
				if(xmlLocation == std::string::npos)
				{
					xmlLocation = tempFileName.find(".XML");
				}

				// Not found, try this one mixed case extension
				if(xmlLocation == std::string::npos)
				{
					xmlLocation = tempFileName.find(".Xml");
				}

				// Finally, we have never found a .xml extension, so add one
				if( xmlLocation == std::string::npos )
				{
					tempFileName+= ".xml";
				}
                os.open(tempFileName.c_str(), ios::out);
                if(os.is_open())
                  {
                    XMLFileName = tempFileName;
                  }
              }
            else
              {
                os.open(XMLFileName.c_str(), ios::out);
              }
          } while( !os.is_open() );

        std::cout << "Outputting to file: " << XMLFileName << std::endl;

          // Output our XML file.  Yes this should use a library, don't shoot the maintenance coder.
          if( os.is_open() )
	{	/****************************************************************************************************************************/
		// start printing to the file according the the file format
		os << "<Routing_Trees start_time=\"" << startTime.c_str() << "\" end_time=\"" << endTime.c_str() << "\" units=\"unix_time_ms\">\n\t<Nodes>";		
                // start printing the tree nodes information here
		vector<int>	indexToFetchNodeIndexVec;							// for each node in the tree, it will have an index value used for outputing
		int treeNodeIndex = 0;											// record all the tree nodes one by one starting from 1
		for(int i=0; i<entries.size()+nodes.size()+fixes.size(); i++)
		{
			Node* temp = fetchNode(i);
			if(temp->treeNodeOrNot()) // if we have a tree node
			{
				// the position in indexToFetchNodeIndexVec is the index value of the node, the stored value is the one used in fetchNode(int i) function
				indexToFetchNodeIndexVec.push_back(i);
				treeNodeIndex++;
				os << "\n\t\t<Node index=\"" << treeNodeIndex << "\" latitude=\"" << centerLati+temp->getX()*latiPerPixel << "\" longitude=\""
				 << centerLong+temp->getY()*longPerPixel << "\" altitude=\"" << ALTITUDE_AT_BASE_PLANE+ALTITUDE_PER_PIXEL*temp->getZ() << "\" alt_unit=\"ft\" ";
				switch(temp->getNodeType())								// print the node type
				{
					case ENTRY_NODE:
						os << "type=\"source\">";
						break;
					case INTERNAL_NODE:
						os << "type=\"internal\">";
						break;
					case FIX_NODE:
						os << "type=\"sink\">";
						break;
				}
				os << "\n\t\t\t<Operational_Flexibility_Discs>";					// print the operational flexity pairs
				for(int j=0; j<temp->getFreeRadiusVecSize(); j++)
				{
					double radius, prob;
					temp->getFreeRadiusResults(j, &radius, &prob);
					os << "\n\t\t\t\t<Disc index=\"" << j+1 << "\" radius=\"" << radius*NMILESPERPIXEL << "\" units=\"nm\" probability=\"" << prob << "\"/>";
				}
				os << "\n\t\t\t</Operational_Flexibility_Discs>\n\t\t</Node>";		
			}
		}
		// output the information for deviation nodes of edges. the deviationNodeIndex variable is for edges' deviation nodes
		int deviationNodeIndex = treeNodeIndex;						
		// output all the operational flexibility deviation nodes information for each edge
		for(int i=0; i<edges.size(); i++)
		{
			if(edges[i]->treeEdgeOrNot())
			{
				for(int j = 0; j<edges[i]->getDeviationNodesSize(); j++)		// output each deviation node information
				{
					Node* temp = edges[i]->getDeviationNode(j);
					deviationNodeIndex++;
					os << "\n\t\t<Node index=\"" << deviationNodeIndex << "\" latitude=\"" << centerLati+temp->getX()*latiPerPixel;
					os << "\" longitude=\"" << centerLong+temp->getY()*longPerPixel << "\" altitude=\"" << ALTITUDE_AT_BASE_PLANE+ALTITUDE_PER_PIXEL*temp->getZ();
					os << "\" alt_unit=\"ft\" ";
					os << "type=\"deviation\">";
					os << "\n\t\t\t<Operational_Flexibility_Discs>";		// print the operational flexity pairs
					for(int j=0; j<temp->getFreeRadiusVecSize(); j++)
					{
						double radius, prob;
						temp->getFreeRadiusResults(j, &radius, &prob);
						os << "\n\t\t\t\t<Disc index=\"" << j+1 << "\" radius=\"" << radius*NMILESPERPIXEL << "\" units=\"nm\" probability=\"" << prob << "\"/>";
					}
					os << "\n\t\t\t</Operational_Flexibility_Discs>\n\t\t</Node>";	
				}
			}
		}
		os << "\n\t</Nodes>";
		// now after the nodes, outputting the edges information
		vector<int>	indexToTreeEdgeIndexVec;								// for each edge in the tree, it will have an index value used for outputing
		int treeEdgeIndex = 0;												// record all the tree nodes one by one starting from 1
		deviationNodeIndex = treeNodeIndex;									// the last output index of tree nodes
		os << "\n\t<Arcs>";
		for(int i=0; i<edges.size(); i++)
		{
			if(edges[i]->treeEdgeOrNot())
			{
				// get the head node's printed index, in the indexToFetchNodeIndexVec vector, where the position in the vector is the index-1
				int unifiedIndexHead = getFetchNodeIndex(edges[i]->getHead());	// the head's unified index
				int tempIndexHead = find(indexToFetchNodeIndexVec.begin(), indexToFetchNodeIndexVec.end(), unifiedIndexHead)-indexToFetchNodeIndexVec.begin()+1;
				int unifiedIndexTail = getFetchNodeIndex(edges[i]->getTail());	// the head's unified index
				int tempIndexTail = find(indexToFetchNodeIndexVec.begin(), indexToFetchNodeIndexVec.end(), unifiedIndexTail)-indexToFetchNodeIndexVec.begin()+1;
				indexToTreeEdgeIndexVec.push_back(i);						// the position in the edges vector to its printed index
				treeEdgeIndex++;
				os << "\n\t\t<Arc index=\"1" << treeEdgeIndex << "\" start_node=\"" << tempIndexHead << "\" end_node=\"" << tempIndexTail << "\" cost=\"";
				os << edges[i]->getLength()*NMILESPERPIXEL << "\" cost_unit=\"nm\">";
				// print rnp operational flexity pairs
				os << "\n\t\t\t<RNP_Levels>";
				for(int j=0; j<edges[i]->getRNPVecSize(); j++)
				{
					double radius, prob;
					edges[i]->getRNPResults(j, &radius, &prob);
					os << "\n\t\t\t\t<RNP_Level RNP=\"" << radius*NMILESPERPIXEL << "\" probability=\"" << prob << "\"/>";
				}
				os << "\n\t\t\t</RNP_Levels>";
				// print path streching information on the right side of the edge
				for(int j=0; j<edges[i]->getPathStretchingVecSize(); j++)
				{
					double radius, prob;
					edges[i]->getPathStretchingResults(j, &radius, &prob);
					os << "\n\t\t\t<Operational_Flexibility_Rectangle index=\"" << j+1 << "\" width=\"" << radius*NMILESPERPIXEL << "\" units=\"nm\" ";
					os << "probability=\"" << prob << "\""; 
					os << " position=\"right\" type=\"path_stretch\"/>\n\t\t\t<Rectangle_Coord direction=\"clockwise\">";
					/****************************************************************/
					// print the rectangle information here in clock wise order
					double x1, y1, x2, y2;									// compute the vertices on the right of the current edge
					edges[i]->computeRightSideRectangleVertices(radius, &x1, &y1, &x2, &y2);
					os << "\n\t\t\t\t\t<Coord index=\"1\" lat=\"" << centerLati+latiPerPixel*edges[i]->getHead()->getX() << "\" lon=\"";
					os << centerLong+longPerPixel*edges[i]->getHead()->getY() << "\"/>";
					os << "\n\t\t\t\t\t<Coord index=\"2\" lat=\"" << centerLati+latiPerPixel*edges[i]->getTail()->getX() << "\" lon=\"";
					os << centerLong+longPerPixel*edges[i]->getTail()->getY() << "\"/>";
					os << "\n\t\t\t\t\t<Coord index=\"3\" lat=\"" << centerLati+latiPerPixel*x1 << "\" lon=\"" << centerLong+longPerPixel*y1 << "\"/>";
					os << "\n\t\t\t\t\t<Coord index=\"4\" lat=\"" << centerLati+latiPerPixel*x2 << "\" lon=\"" << centerLong+longPerPixel*y2 << "\"/>";
					/****************************************************************/
					os << "\n\t\t\t</Rectangle_Coord>";
				}
				// print wiggle room information on the left side of the edge
				for(int j=0; j<edges[i]->getWiggleRoomVecSize(); j++)
				{
					double radius, prob;
					edges[i]->getWiggleRoomResults(j, &radius, &prob);
					os << "\n\t\t\t<Operational_Flexibility_Rectangle index=\"" << j+1 << "\" width=\"" << radius*NMILESPERPIXEL << "\" units=\"nm\" probability=\"" << prob << "\""; 
					os << " position=\"left\" type=\"holding\"/>\n\t\t\t<Rectangle_Coord direction=\"clockwise\">";
					/****************************************************************/
					// print the rectangle information here in clock wise order
					double x1, y1, x2, y2;									// compute the vertices on the left of the current edge
					edges[i]->computeLeftSideRectangleVertices(radius, &x1, &y1, &x2, &y2);
					os << "\n\t\t\t\t\t<Coord index=\"1\" lat=\"" << centerLati+latiPerPixel*x1 << "\" lon=\"" << centerLong+longPerPixel*y1 << "\"/>";
					os << "\n\t\t\t\t\t<Coord index=\"2\" lat=\"" << centerLati+latiPerPixel*x2 << "\" lon=\"" << centerLong+longPerPixel*y2 << "\"/>";
					os << "\n\t\t\t\t\t<Coord index=\"3\" lat=\"" << centerLati+latiPerPixel*edges[i]->getTail()->getX() << "\" lon=\"";
					os << centerLong+longPerPixel*edges[i]->getTail()->getY() << "\"/>";
					os << "\n\t\t\t\t\t<Coord index=\"4\" lat=\"" << centerLati+latiPerPixel*edges[i]->getHead()->getX() << "\" lon=\"";
					os << centerLong+longPerPixel*edges[i]->getHead()->getY() << "\"/>";
					/****************************************************************/
					os << "\n\t\t\t</Rectangle_Coord>";
				}
				// print the deviation nodes along an edge
				os << "\n\t\t\t<Off_Nominal_Exit_Points>";
				for(int j=0; j<edges[i]->getDeviationNodesSize(); j++)
				{
					deviationNodeIndex++;
					os << "\n\t\t\t\t<Off_Nominal_Exit_Point index=\"" << deviationNodeIndex << "\"/>";
				}
				os << "\n\t\t\t</Off_Nominal_Exit_Points>\n\t\t</Arc>";
			}
		}
		os << "\n\t</Arcs>";
		// output tree branch information here. For each tree branch, output its edge lists
		os << "\n\t<Branches>";
		for(int i=0; i<entries.size(); i++)
		{
			if(entries[i]->getDrawingRNP()==0)					// make sure that the branch has actual demand
				continue;
			Node* temp = entries[i];
			os << "\n\t\t<Branch index=\"" << 100+i << "\" rnp_equipage=\"" << temp->getDrawingRNP()*NMILESPERPIXEL << "\">";
			while(temp->getNodeType()!=FIX_NODE)
			{
				Edge* tempEdge = temp->getOutEdge(temp->getTreeOutEdgeIndex());
				// the original index of the edge in the original edge vector
				int edgeOriginalIndex = find(edges.begin(), edges.end(), tempEdge)-edges.begin();
				int printingIndex = find(indexToTreeEdgeIndexVec.begin(), indexToTreeEdgeIndexVec.end(), edgeOriginalIndex)-indexToTreeEdgeIndexVec.begin()+1;
				os << "\n\t\t\t<Branch_arc index=\"1" << printingIndex << "\"/>";
				temp = temp->getOutNode(temp->getTreeOutEdgeIndex());
			}
			os << "\n\t\t</Branch>";
		}
		os << "\n\t</Branches>\n</Routing_Trees>\n";
		// finish printing, close the file stream and return file generated successfully
		/****************************************************************************************************************************/
		os.close();
		return true;
	}
	std::cerr << "\nCannot successfully create the file..." << std::endl;
	return false;
}

// find a node based on its layer and layerIndex, return NULL means no matching node is found
Node* RoutingDAG::findNode(int layer, int layerIndex)
{
	// search all the nodes linearly to find the correnponding one
	for(int i=0; i<nodes.size(); i++)
		if(nodes[i]->getLayer()==layer && nodes[i]->getLayerIndex()==layerIndex)
			return nodes[i];
	for(int i=0; i<fixes.size(); i++)
		if(fixes[i]->getLayer()==layer && fixes[i]->getLayerIndex()==layerIndex)
			return fixes[i];
	for(int i=0; i<entries.size(); i++)
		if(entries[i]->getLayer()==layer && entries[i]->getLayerIndex()==layerIndex)
			return entries[i];
	return NULL;
}

// after all the nodes are read in, generate the edges for the nodes based on the angle requirements
// centerAngle is a standard value that the edges must have angles with in 60 degrees(PI/3) of the center angle
bool RoutingDAG::generateEdgeSet()
{
	if(nodesReadIn == NODES_NOT_READ_IN)
	{
		std::cerr << "\nThe Nodes have to be read in before generating the edges." << std::endl;
		return false;
	}
	/***************************************************************************************************/
	// generate the edge set base on the nodes that are read in. Each node has outgoing edges going into the nodes in the next 3 layers, as long as the edge
	// is with in PI/3 degrees of the centerAngle, the order of the edges are from rightmost to leftmost, from furthest to nearest, except for the nodes that have
	// edges going to fix nodes, then the order would be from leftmost to rightmost
	if(!generateLayerStartingIndexVector())					// generate the layer starting position vector first
		return false;
	// iterate for each layer, deal with the fix node layer and the second to the last layer separately later because the last layer only has edges going to fix nodes
	// but the order where we fill in fix nodes is opposite
	for(int i=0; i<numLayers-2; i++)					
	{
		int startingLayer = (i+3>numLayers-2)? (numLayers-2) : (i+3);
		int endingLayer = i+1;							// connect to points that lie in the range [startingLayer, endingLayer]
		// loop thru all the nodes int the current layer
		for(int j = layerStartingIndex[i]; j<layerStartingIndex[i+1]; j++)		
		{
			Node* current = fetchNode(j);				// the current Node on the current Layer
			double x = current->getX();
			double y = current->getY();
			int *shortestDistanceNodeIndex = new int[startingLayer-endingLayer+1];
			// compute for each layer(maximum 3), the closest node on that layer to the current node
			for(int k=endingLayer; k<=startingLayer; k++)
			{
				// the distance from current node to the first node on the kth layer
				Node* temp = fetchNode(layerStartingIndex[k+1]-1);
				double shortestDistance = (x-temp->getX())*(x-temp->getX()) + (y-temp->getY())*(y-temp->getY());
				shortestDistanceNodeIndex[k-endingLayer] = layerStartingIndex[k+1]-1;
				for(int l = layerStartingIndex[k+1]-2; l>=layerStartingIndex[k]; l--)
				{
					temp = fetchNode(l);
					double tempDist = (x-temp->getX())*(x-temp->getX()) + (y-temp->getY())*(y-temp->getY());
					if(tempDist < shortestDistance)		// compute the closest distance and record its position
					{
						shortestDistance = tempDist;
						shortestDistanceNodeIndex[k-endingLayer] = l;
					}
				}
			}			
			// in the shortestDistanceNodeIndex array, we now have the closest node on that specific layer to the current node
			// we then generate edges, from endingLayer to startingLayer, from "right" to "left", at most 15 edges(5 to each layer) are generated for each node
			for(int k=2; k>=-2; k--)					// at most 5 nodes on each levels 
			{
				for(int l = startingLayer; l>=endingLayer; l--)
				{
					if(l==endingLayer && (k==2 || k==-2))									// the next layer has only 3 nodes
						continue;
					int tempIndex = shortestDistanceNodeIndex[l-endingLayer]+k;				// we are interested in this node at this step
					if(tempIndex>=layerStartingIndex[l+1] || tempIndex<layerStartingIndex[l])		// a node in the next layer or previous layer, then not feasible
						continue;
					Node* tempNode = fetchNode(tempIndex);									// this node will form one of the edges with the current node
					Edge* tempEdge = new Edge(current, tempNode);
					current->insertOutNodeEdge(tempNode, tempEdge);							// insert the node and edge into the neighbor vector of the current Node
					tempNode->insertInNodeEdge(current, tempEdge);							// the other node points back to the current node, too
					edges.push_back(tempEdge);
				}
			}
			delete []shortestDistanceNodeIndex;
		}
	}
	/***************************************************************************************************/
	// deal with the connection to the final layer: the layer that holds fix nodes
	// we suppose that only the last 2 layers connect to the fix nodes, so we add the fix nodes into their outgoing edge list one by one
	int startingLayer = numLayers-3>=0? numLayers-3 : 0;
	for(int i=layerStartingIndex[startingLayer]; i<layerStartingIndex[numLayers-1]; i++)		// these nodes will be connected to all fix nodes
	{
		Node* current = fetchNode(i);
		int fixNodesSize = fixes.size();
		for(int j=0; j<fixNodesSize; j++)
		{
			Node* tempNode = fixes[j];														// connect current node to the fix node, from small angle to large angle
			Edge* tempEdge = new Edge(current, tempNode);									// in the opposite order compared to previous inserting order
			current->insertOutNodeEdge(tempNode, tempEdge);
			tempNode->insertInNodeEdge(current, tempEdge);	
			edges.push_back(tempEdge);
		}
	}
	/***************************************************************************************************/
	edgesGenerated = EDGES_GENERATED;
	return true;
}

// a helper function to fill in the layerStartingIndex vector, which represented the starting position of each layer of nodes
bool RoutingDAG::generateLayerStartingIndexVector()
{
	if(nodesReadIn == NODES_NOT_READ_IN)
	{
		std::cerr << "\nThe Nodes have to be read in before generating the edges." << std::endl;
		return false;
	}
	layerStartingIndex.push_back(0);						// the starting position of the entry nodes is firstly pushed in the vector
	int initialLayer = 0;
	int entrySize = entries.size();
	for(int i=0; i<nodes.size(); i++)
	{
		if(nodes[i]->getLayer()!=initialLayer)
		{
			layerStartingIndex.push_back(i+entrySize);
			initialLayer = nodes[i]->getLayer();
		}		// now the layerStartingIndex vector stores the starting position of each layer in the nodes vector
	}
	layerStartingIndex.push_back(entrySize+nodes.size());	// the starting position of the fix nodes is lastly pushed in the vector
	layerStartingIndex.push_back(entrySize+nodes.size()+fixes.size());		// the last position, after the end of fix nodes
	return true;
}

// for a node, get its index value used in the fetchNode function as parameter(kind of a unified index among all the nodes)
int RoutingDAG::getFetchNodeIndex(Node* temp)
{
	if(temp->getNodeType()==ENTRY_NODE)
	{
		if(find(entries.begin(), entries.end(), temp)!=entries.end())		// if the node does exist
			return find(entries.begin(), entries.end(), temp)-entries.begin();
	}
	else if(temp->getNodeType()==INTERNAL_NODE)
	{
		if(find(nodes.begin(), nodes.end(), temp)!=nodes.end())				// if the node does exist
			return find(nodes.begin(), nodes.end(), temp)-nodes.begin()+entries.size();
	}
	else if(temp->getNodeType()==FIX_NODE)
	{
		if(find(fixes.begin(), fixes.end(), temp)!=fixes.end())				// if the node does exist
			return find(fixes.begin(), fixes.end(), temp)-fixes.begin()+entries.size()+nodes.size();
	}
	return -1;
}


// find a node based on its index among the big vector entries, nodes and fixes
Node* RoutingDAG::fetchNode(int n)
{
	int entrySize = entries.size();
	int nodesSize = nodes.size();
	int fixSize = fixes.size();
	if(n>=0 && n<entrySize)
		return entries[n];
	if(n>=entrySize && n<nodesSize+entrySize)
		return nodes[n-entrySize];
	if(n>=nodesSize+entrySize && n<nodesSize+entrySize+fixSize)
		return fixes[n-entrySize-nodesSize];
	return NULL;
}

/***************************************************************************************************************************************************************/
// the key function, routing algorithm, used to generate a tree from a prebuilt graph, given the weather, rnp for each entry point and threshold value
// generating a tree by marking the tree edges and nodes their flag values to mark them as in the tree
// effectiveThres means if a weather cell's deviation probability is below this value, it's going to be considered as NULL
// routingThres means that we compute the weighted total probability of the weathercells p1*(0 or 1) + p2*(0 or 1) +..., 
// if the value < routingThres, then it is considered an obstacle 
bool RoutingDAG::generateTree(const vector<WeatherData> &wData, vector<float> rnp, float effectiveThres, float routingThres)
{
	// before computing, first reset all tree related variables and tree generating status
	resetTree();	
	status = TREE_NOT_GENERATED;	
	treeShapeStatus = BOTTOM_TREE;
	if(nodesReadIn == NODES_NOT_READ_IN || edgesGenerated == EDGES_NOT_GENERATED)
	{
		std::cerr << "\nNodes or Edges are not ready yet to generate a tree." << std::endl;
		return false;
	}
	if(rnp.size()!=entries.size())							// input error, there is supposed to be an rnp value related to each entry point
	{														// format error, then prompt the problem
		std::cerr << "\nNumber of RNP values does not match the number of entry points." << std::endl;
		return false;
	}
	for(int i=0; i<entries.size(); i++)						// first make sure that all the entry nodes are themselves weather obstacle free
		if(entries[i]->testRadiusWithWeatherDataSet(rnp[i], wData, effectiveThres, routingThres))
		{
			return false;	
		}
	// then start testing branches coming out of each entry node
	for(int i=0; i<entries.size(); i++)
	{
		// test if the entry node is feasible, then route a branch for each entry node
		for(int j=0; j<entries.size()+nodes.size()+fixes.size(); j++)
		{
			fetchNode(j)->setVisited(NOT_VISITED);					// when starting a new graph search from a new entry node, set every node unvisited
		}
		// first update the layerUsedIndex vector by walking thru the previous branches
		updateLayerUsedIndexVector(i);
		entries[i]->setDrawingRNP(rnp[i]);							// set the rnp information and treeNode status of the entry node
		entries[i]->setTreeNode();
		for(int j=0; j<entries.size()+nodes.size()+fixes.size(); j++)
		{
			fetchNode(j)->setVisited(NOT_VISITED);					// when starting a new search, all the nodes within the search DAG are set to be unvisited
		}
		if(!routeBranch(entries[i], i, wData, rnp[i], effectiveThres, routingThres))		
		{
			resetTree();											// tree routing failed, then no tree is generated, reset
			return false;
		}
		layerUsedIndex[0] = i;										// the current node is routed the entry level, which nodes were used already
		// process the tree brach from the current entry node so that each tree node in the branch knows which its outgoing edge in the tree is
		// the function sets the treeOutEdgeIndex variable
		treeBranchPostProcessing(entries[i], rnp[i]);
	}
	status = TREE_GENERATED;
	treeShapeStatus = BOTTOM_TREE;
	return true;
}

// search from a specific entry node "start", using DFS to find a path to a fix node
// effectiveThres means if a weather cell's deviation probability is below this value, it's going to be considered as NULL
// routingThres means that we compute the weighted total probability of the weathercells p1*(0 or 1) + p2*(0 or 1) +..., 
// if the value < routingThres, then it is considered an obstacle 
bool RoutingDAG::routeBranch(Node *start, int entryIndex, const vector<WeatherData> &wData, float rnp, float effectiveThres, float routingThres)
{
	// if there is no entry demand at all from this entry node, then simply return true
	if(rnp==0)	return true;
	/**************************************************************************************************************************/
	stack<Node*> tempStack;											// a stack that is used to do the DFS starting at node start
	tempStack.push(start);											// push the start node itself into the stack
	// do the loop while the stack is not empty, doing a Depth First Search
	while(!tempStack.empty())
	{
          std::cout  <<  "."  <<  std::flush;		
		Node* temp = tempStack.top();								// pop the first node out and store it as Node* temp
		tempStack.pop();
		for(int i=0; i<=temp->getOutSize()-1; i++)					// test each of its outgoing edges (from top(right) to bottommost(left))
		{															// because the stack pops elements last in first out
			Edge* tempEdge = temp->getOutEdge(i);					// the pointers to the current outer edge and outer node
			Node* tempNode = temp->getOutNode(i);
			// if the new outgoing node/edge pair conflict with weather, then just ignore
			if(tempEdge->testRNPWithWeatherDataSet(rnp, wData, effectiveThres, routingThres) ||
			   tempNode->testRadiusWithWeatherDataSet(rnp, wData, effectiveThres, routingThres))
				continue;
			// we explore the node only if it's on the right side of the existing tree(bottommost filling, so cannot cross previous branches)
			if(tempNode->getLayerIndex()>layerUsedIndex[tempNode->getLayer()])
			{
				// if the node or edge conflicts with the previous branch, then it must not be used and means all the nodes "below" the node are not useful anymore
				if(testPreviousBranchNode(tempNode, rnp, entryIndex) || testPreviousBranchEdge(tempEdge, rnp, entryIndex))				
					continue;
				if(tempNode->getInDegree()<2)						// if we find a node that is feasible, and has indegree less than 2
				{
					tempNode->setPrevTreeNode(temp);						// mark the prevNode and prevEdge in the tree structure
					tempNode->setPrevTreeEdge(tempEdge);
					// now we push this node into the stack, if the node was never visited before, then push it directly into the stack
					if(tempNode->getVisited()==NOT_VISITED)
					{
						tempNode->setVisited(VISITED);
						tempStack.push(tempNode);						// if not visited before, push all the feasible nodes into the stack
					}
					// but this node might already be in the stack, so in this case, we pop it out first, then push it into the stack again so that it's now the top element
					else if(tempNode->getVisited()==VISITED)			// test if it's already in the stack
					{
						stack<Node*> tempBufferStack;
						bool nodeInsideStack = false;
						while(!tempStack.empty())
						{
							Node* tempTransferingNode = tempStack.top();
							tempStack.pop();
							if(tempTransferingNode == tempNode)			// if we find the node in the stack, setting the boolean variable and jump out of the loop
							{
								nodeInsideStack = true;
								break;
							}
							else
								tempBufferStack.push(tempTransferingNode);		// temporarily storing the node into another stack
						}
						while(!tempBufferStack.empty())					// move the elements back into the original stack
						{
							tempStack.push(tempBufferStack.top());
							tempBufferStack.pop();
						}
						if(nodeInsideStack)								// means we are going to push this node in the stack again
							tempStack.push(tempNode);
					}
				}
			}
			/**********************************************************************************************************************************************************/
			// THE 2 ENDING CONDITIONS OF FINDING A BRANCH, THESE 2 CASES MARKS THE END OF A BRANCH
			// this case happens then we are merging into an existing branch(but not at the fix node). 
			// We test all the remaining part of the branch to see if the remaining branch is feasible, 
			// and make sure that this node is from the immediate previous branch, instead of some node from other previous branches
			// note that once merged into the brach of the tree, it cannot deviate out from the branch 
			// and this never happens for the first branch of the tree, because layerUsedIndex vector is all set to -1
			if(tempNode->getNodeType()!=FIX_NODE &&	tempNode->getLayerIndex()==layerUsedIndex[tempNode->getLayer()] && 
			   tempNode->treeNodeOrNot() && tempNode->getInDegree()<2 && onPreviousBranch(tempNode, entryIndex) &&
			   !testPreviousBranchTillCurrentLayerEdge(tempEdge, rnp, entryIndex))		// test if the edge collide with part of previous branch (till current layer)
			{																			// test the minimum distance between 2 merging nodes constraint
				if(testDistanceTooCloseToMergingNodesOnPreviousBranch(tempNode, entryIndex))
					continue;
				// meaning that the remaining part of the tree would work, then just finish building the tree
				if(testRemainingBranchWhileMerging(tempNode, wData, rnp, effectiveThres, routingThres))		
				{
					setTreeBranchUpMerging(tempNode, rnp);									// set the drawingRNP values from the fix node to the node that we start
					tempNode->setPrevTreeNode(temp);									// mark the prevNode and prevEdge in the tree structure
					tempNode->setPrevTreeEdge(tempEdge);
					setTreeBranchUp(tempNode, start, rnp);								// set the tree information and the branch is complete, return true, find a branch
					return true;
				}
			}
			// if we reach a fix node, then mark the branch in the tree, meaning that we just found a complete branch
			if(tempNode->getLayerIndex()>=layerUsedIndex[tempNode->getLayer()] && tempNode->getNodeType()==FIX_NODE)			
			{	// if we are approaching a merging fix node here, then test the min distance constraint first
				if(tempNode->getInDegree()==1 && testDistanceTooCloseToMergingNodesOnPreviousBranch(tempNode, entryIndex))	
					continue;
				if(testPreviousBranchTillCurrentLayerEdge(tempEdge, rnp, entryIndex))
					continue;
				// no need to test the node itself again since it's a fix node, already tested when generated
				if(tempNode->getInDegree()<2)
				{
					tempNode->setPrevTreeNode(temp);
					tempNode->setPrevTreeEdge(tempEdge);
					// set all the edges that we used to be in the tree, and it is guaruateed to be a bottommost tree branch
					// trace back in the tree and mark everything as in the tree
					setTreeBranchUp(tempNode, start, rnp);
					// the startNode is already in the tree, therefore no need to loop till that node
					return true;								// means we find a path starting from the current entry node
				}
			}	/**************************************************************************************************************************************************/
		}
	}
	return false;				// after exploring all the possible branches from this particular node, there is nothing feasible, report NO branch found
}

// test a branch of the tree starting from a given Node, and see if the remaining part of the branch is good to be merged in by another branch
bool RoutingDAG::testRemainingBranchWhileMerging(Node *start, const vector<WeatherData> &wData, float rnp, float effectiveThres, float routingThres)
{
	if(start->getInDegree()>=2)										// if a node is already fully occupied, then return fail
		return false;
	Node* temp = start;												// the node used to walk down the remaining part of the tree
	while(temp)
	{
		// getDrawingRNP is used to denote if the previous RNP is alreayd larger, then no need to test again
		// test if the current node is clear of weather obstacle with the new rnp				
		if(temp->getDrawingRNP()<rnp && temp->testRadiusWithWeatherDataSet(rnp, wData, effectiveThres, routingThres))	
			return false;
		if(temp->getNodeType()==FIX_NODE)							// if we are getting to a fix node, then we suceessfully prove that the branch is safe to be merged in
			break;													
		Edge* tempEdge = temp->getOutEdge(temp->getTreeOutEdgeIndex());
		if(tempEdge)		// if we find a tree edge (each node has at MOST one outgoing edge that is a tree edge), and its NOT clear of weather with the new rnp value
		{
			if(tempEdge->getDrawingRNP()<rnp && tempEdge->testRNPWithWeatherDataSet(rnp, wData, effectiveThres, routingThres))
				return false;									// NOT clear of weather, report false
		}
		else				// means no outgoing tree edge is found, means there is tree error in the previous branch
		{
			std::cerr << "\nWhen merging, detected that the previous branch of the tree is incomplete" << std::endl;
			resetTree();
			return false;
		}
		temp = temp->getOutNode(temp->getTreeOutEdgeIndex());						// test the node on the next level of the tree
	}
	return true;													
}

// set all the edges that we used to be in the tree, and it is guaruateed to be a bottommost tree branch
// when we find the last node of a tree(a fix node), mark the tree from the current node all way up to the entry node start
// set rnp values and in degrees
void RoutingDAG::setTreeBranchUp(Node* current, Node* start, float rnp)
{
	while(current!=start)
	{
		current->addInDegree();				// one more edge is coming in	
		current->setTreeNode();				// the node is a tree node
		if(rnp>current->getDrawingRNP())
			current->setDrawingRNP(rnp);		// the max rnp is updated
		current->getPrevTreeEdge()->setTreeEdge();	// the edge leading to this node is a tree edge
		if(rnp>current->getPrevTreeEdge()->getDrawingRNP())
			current->getPrevTreeEdge()->setDrawingRNP(rnp);
		layerUsedIndex[current->getLayer()] = current->getLayerIndex();		// in current layer, the layerIndex node was the last node used in the tree
		current = current->getPrevTreeNode();
	}
}

// set the drawingRNP values from the fix node to the merging node, only rnp values are updated because the other status variables are already set to be in the tree
// Node* current is always a fix node, and it links up to start node, this funtion is only called when we are merging a new branch into an existing branch
void RoutingDAG::setTreeBranchUpMerging(Node *start, float rnp)
{
	Node* temp = start;								// find the fix node of the branch containing node start
	while(temp->getNodeType()!=FIX_NODE)
		temp = temp->getOutNode(temp->getTreeOutEdgeIndex());
	Node* current = temp;							// current is the fix node of the branch constaining start
	while(current!=start)
	{
		if(rnp>current->getDrawingRNP())
			current->setDrawingRNP(rnp);		// the max rnp is updated
		if(rnp>current->getPrevTreeEdge()->getDrawingRNP())
			current->getPrevTreeEdge()->setDrawingRNP(rnp);
		current = current->getPrevTreeNode();
	}
}

// process a tree brach so that each tree node knows which its outgoing edge in the tree is, by marking its treeOutEdgeIndex variable
void RoutingDAG::treeBranchPostProcessing(Node* start, float rnp)			// the node pointer start always points to an entry node
{
	if(start->getDrawingRNP()==0)								// nothing to process here
		return;
	Node* temp = start;
	while(temp->getNodeType()!=FIX_NODE)
	{
		for(int i=0; i<temp->getOutSize(); i++)
		{
			if(temp->getOutEdge(i)->treeEdgeOrNot())			// if we find a tree edge
			{
				temp->setTreeOutEdgeIndex(i);
				break;			// break out of the for loop so that we get to the next level in the tree
			}
		}
		temp = temp->getOutNode(temp->getTreeOutEdgeIndex());	// next level of tree node
	}			// after the function, if we have the start node, we could use the treeOutEdgeIndex variable to trace the entire tree branch
}

// update the layerUsedIndex vector by walking the previous (from 0 to (i-1)th) branches
void RoutingDAG::updateLayerUsedIndexVector(int entryIndex)
{
	if(entryIndex>=entries.size())					// i should be at most entries.size()-1
		return;
	for(int j=0; j<numLayers; j++)			// first restore every element to -1
		layerUsedIndex[j] = -1;
	for(int i=0; i<entryIndex; i++)
	{
		Node* start = entries[i];
		if(start->getDrawingRNP()==0)		// if there is no branch starting from this entry node
			continue;
		while(start->getNodeType()!=FIX_NODE)
		{
			int tempLayer = start->getLayer();
			layerUsedIndex[tempLayer] = max(layerUsedIndex[tempLayer], start->getLayerIndex());
			start = start->getOutNode(start->getTreeOutEdgeIndex());			// go to the next level of the bottommost tree
		}
		// now start is a fix node
		layerUsedIndex[start->getLayer()] = max(layerUsedIndex[start->getLayer()], start->getLayerIndex());
	}
}

/***************************************************************************************************************************************************************/
// after the bottommost tree is generated, tauten its branches so that it looks much better using Dijkstra algorithm on DAG
// in the front of this function, no feasibility needs to be tested (other than if the bottommost tree was generated already) 
// because the function is only called after the bottommost tree is generated
bool RoutingDAG::generateTautenedTree(const vector<WeatherData> &wData, vector<float> rnp, float effectiveThres, float routingThres)
{
	if(status == TREE_NOT_GENERATED)		// this is not actually going to happen because the same test was always conducted before calling this function
	{
		std::cerr << "\nPlease generate the bottommost tree before tautening its branches!" << std::endl;
		return false;
	}
	else if(treeShapeStatus == TAUTENED_TREE)	// if the tree was generated already, then just return, no need to generate again
		return false;
	// save the original bottommost tree information, in case there are multiple rounds of tautened tree generating
	for(int i=0; i<entries.size()+nodes.size()+fixes.size(); i++)
	{
		fetchNode(i)->storeTreeOutEdgeIndexInformation();
	}
	// when topMostTendency is 1, nothing unusual happens. Then next round, every 2nd branches will try to merge into its next branch
	// next round, then every 3rd and 2nd branches will try to merge into their corresponding 1st branch. until topMostTendency is entry.size(), where
	// each branch will try to merge into the first branch. This is used to deal with the bottleneck cases and urge branches to merge early
	// the tree will look more like a topmost tree if every branch tries to merge into the first branch. Hence it's called topMostTendency
	for(int topMostTendency=1; topMostTendency<=entries.size(); topMostTendency++)
	{
		std::cout  <<  ".";
		//preparation for the new tree: clear the unuseful(useful for the new tree) data structure related to the previous tree for nodes and edges
		for(int i=0; i<entries.size()+nodes.size()+fixes.size(); i++)
		{
			Node* temp = fetchNode(i);
			float tempRNP = temp->getDrawingRNP();
			// set nodes NOT in the tree, in degree 0, BUT THE TREE OUTGOING EDGE INFORMATION IS KEPT(IMPORTANT!!!): we can still trace a branch from an entry
			temp->resetTreeStatus();			
			temp->setDrawingRNP(tempRNP);					// restore the drawingRNP information
			temp->restoreTreeOutEdgeIndexInformation();		// whenever a new round of tautening starts, restore the bottommost tree info
		}
		for(int i=0; i<edges.size(); i++)
			edges[i]->setNotTreeEdge();						// all edges are not NOT in the tree
		// note that we can still trace a branch from any entry nodes using the treeOutEdgeIndex variable
		//set these indices one more than the index of the last node in the layer, so no node is used at the beginning
		layerUsedIndexReverseDirection.clear();
		for(int i=0; i<numLayers; i++)
			layerUsedIndexReverseDirection.push_back(layerStartingIndex[i+1]-layerStartingIndex[i]);
		// routing from the entry nodes, from last one to the first one
		bool successfullyTautened = true;					// denote if we route successfully till the last entry node
		for(int i=entries.size()-1; i>=0; i--)
		{
			for(int j=0; j<entries.size()+nodes.size()+fixes.size(); j++)
				fetchNode(j)->setVisited(NOT_VISITED);		// when starting a new search, all the nodes within the search DAG are set to be unvisited
			for(int j=0; j<entries.size()+nodes.size()+fixes.size(); j++)
				fetchNode(j)->setDistance(*max_element(rnp.begin(), rnp.end())*5*numLayers);	// at the beginning of Dijkstra, every node has very large distance
			// first update the layerUsedIndex vector by walking thru the previous branches
			updateLayerUsedIndexVector(i);
			entries[i]->setTreeNode();
			// if tautening fails, reason is that maybe there is bottleneck place that branch (i+1) takes too much freedom of branch (i), such as
			// a merging node is taken, or the min distance between merging nodes requirement cannot be fulfilled, etc...
			if(!routeTautenedTreeBranch(entries[i], i, wData, rnp[i], effectiveThres, routingThres, topMostTendency))	
			{
				if(topMostTendency==entries.size())			// if the last try also failed, then there is no tautening
				{
					std::cerr << "Tautening algorithm Error! Cannot tauten the tree!" << std::endl;
					resetTree();	
					return false;
				}
				successfullyTautened = false;				// go directly to the next round of tautening with topMostTendency++
				break;
			}
			else								// the tree was generated successfully, then set the corresponding status variables 
				treeBranchPostProcessingForTreeTautening(entries[i], rnp[i]);
		}
		if(successfullyTautened)							// if this round of topMostTendency gives us a tautened tree, then stop
			break;
	}
	// the tree was generated successfully(otherwise returned in the if block), then set the corresponding status variables 
	status = TREE_GENERATED;
	treeShapeStatus = TAUTENED_TREE;
	return true;
}

// the function is only used to compare 2 nodes when stl::sort is used on the deque data structure for shortest path algorithm
bool compareNodes(Node* n1, Node* n2)
{
	if(n1->getLayer()<n2->getLayer())
		return true;
	return false;
}

// this function route a branch of the tautened tree starting from an entry node, the topMostTendency parameter denotes if the tree will more like a 
// topmost tree or not(sometimes the only way to do this is a topmost tree), topMostTendency is from 1 to entries.size() 
bool RoutingDAG::routeTautenedTreeBranch(Node *start, int entryIndex, const vector<WeatherData> &wData, float rnp, 
										 float effectiveThres, float routingThres, int topMostTendency)
{
	if(start->getDrawingRNP()==0)			// if there is no demand from the current entry node, then no branch is coming out of it, return
		return true;
	deque<Node*> tempQueue;					// the queue used for breath first search in Dijkstra algorithm for a DAG
	// decide which is the node to end the current branch ( we may meet at most 3 fix nodes, the one we pick is the one that has the minimum distance
	Node* endingFixNode = NULL;													
	tempQueue.push_back(start);		
	start->setDistance(0);					// the source vertex has distance 0
	// the pointer to the node that the current branch can first merge into the next branch. See more at the last if statement in the while loop
	Node* earlyMergeNode = NULL;			
	while(!tempQueue.empty())
	{
		sort(tempQueue.begin(), tempQueue.end(), compareNodes);		// sort the queue so that the nodes with lower layer comes out first. DAG shortest path algorithm
		Node* temp = tempQueue.front();				// pop the first node of the queue and point to it by temp Pointer
		tempQueue.pop_front();
		/****************************************************************************************************************/
		/* this happens only when we already have a merging node to merge into the "next" branch. We end the current branch routing here because
		 when temp is earlyMergeNode, we are sure that we have the minimum distance from the current entry node to earlyMergeNode.
		 More details see the last if statement in the while loop*/
		if(earlyMergeNode && temp==earlyMergeNode)
		{
			setTreeBranchUpMergingForTreeTautening(earlyMergeNode, rnp);
			setTreeBranchUpForTreeTautening(earlyMergeNode, start, rnp);
			return true;
		}
		/****************************************************************************************************************/
		for(int i=0; i<temp->getOutSize(); i++)		// iterate thru each adjacent node and edge of the current node
		{
			Node* tempNode = temp->getOutNode(i);
			Edge* tempEdge = temp->getOutEdge(i);
			// first test if the new pair of nodes and edges collide with the weather data, if so, then just ignore this pair
			if(tempNode->testRadiusWithWeatherDataSet(rnp, wData, effectiveThres, routingThres) ||
			   tempEdge->testRNPWithWeatherDataSet(rnp, wData, effectiveThres, routingThres))
				continue;							// the node is infeasible, look at the next node
			// then test if this new edge obviously goes across the previous or next branch, then ignore as well
			if(tempNode->getLayerIndex()>layerUsedIndexReverseDirection[tempNode->getLayer()] || 
			   tempNode->getLayerIndex()<layerUsedIndex[tempNode->getLayer()])
				continue;
			// if the edge was NOT completely a part of the next branch (if it is, we don't care about the indegree) 
			if(!(onNextBranch(temp, entryIndex) && onNextBranch(tempNode, entryIndex))		// not going along the next branch
			   && tempNode->getInDegree()==2)		// if the node is currently indegree full, then it's not going to be used again, ignore...	
				continue;
			if(onNextBranch(tempNode, entryIndex) && !onNextBranch(temp, entryIndex))		// if merging into the next branch
			{
				if(testnextBranchTillCurrentLayerEdge(tempEdge, rnp, entryIndex))			// if it crosses part of the next branch(with larger entry index)
					continue;
				if(testDistanceTooCloseToMergingNodesOnNextBranch(tempNode, entryIndex))	// if it does NOT satisfy the min distane between merging nodes constraint
					continue;
				// if the new rnp is too large when merging, then we cannot used this node as a merging node(once merged, cannot depart later)
				if(!testRemainingBranchWhileMerging(tempNode, wData, rnp, effectiveThres, routingThres))	
					continue;
			}
			// if the start of the current edge is already a part of the previous branch, then the new node must also be a node in the previous branch
			// the reason is that once two branches merges together, they have to be together till the fix node
			if(onNextBranch(temp, entryIndex) && tempNode!=nextNodeOfGivenNodeOnNextBranch(temp, entryIndex))
				continue;
			// the new tree branch sits in between the next branch and the previous branch in the bottommost tree, it is allowed to use nodes
			// in the previous and next branches as long as they are not full in indegree(for next branch only)
			bool pushIntoQueue = false;				// if the node is finally going to be pushed into the queue
			// the next are the 3 cases that we consider the node tempNode is a feasible addition to the tree nodes
			if(onNextBranch(tempNode, entryIndex))	// if the node is on next branch (has layerIndex number larger than the current)
				pushIntoQueue = true;
			// if the node is NOT on the next or previous branch, and does not collide with the next branch or previous branch
			else if(!onNextBranch(tempNode, entryIndex) && !onPreviousBranch(tempNode, entryIndex) && !testNextBranchNode(tempNode, rnp, entryIndex) && 
					!testNextBranchEdge(tempEdge, rnp, entryIndex) && !testPreviousBranchNode(tempNode, rnp, entryIndex) && 
					!testPreviousBranchEdge(tempEdge, rnp, entryIndex))	
				pushIntoQueue = true;
			// the new node is allowed to be on the previous branch, but cannot cross it
			else if(onPreviousBranch(tempNode, entryIndex) && !testNextBranchNode(tempNode, rnp, entryIndex) &&!testNextBranchEdge(tempEdge, rnp, entryIndex))
				pushIntoQueue= true;				// (has layerIndex number smaller than the current)
			// if the node is proved to be a feasible choice, push the node into the queue, and update its distance info if the new distance is better
			if(pushIntoQueue)
			{
				if(tempNode->getDistance()>temp->getDistance()+tempEdge->getLength())	// update the distance info
				{
					tempNode->setDistance(temp->getDistance()+tempEdge->getLength());
					tempNode->setPrevTreeNode(temp);
					tempNode->setPrevTreeEdge(tempEdge);
				}
				// note after a node is popped out from the queue, it won't have the chance to go into the queue again, because it's a DAG
				if(tempNode->getVisited()==NOT_VISITED)			// if it's visited before, then only update its distance information
				{												// otherwise, set it visited, and update its distance information
					tempNode->setVisited(VISITED);
					tempQueue.push_back(tempNode);
				}
				// if the node is a feasible FIX node and is closer to the source, then store it using pointer endingFixNode 
				if(tempNode->getNodeType() == FIX_NODE && endingFixNode!=tempNode)
					if (!endingFixNode || (endingFixNode->getDistance()>=tempNode->getDistance()))
						endingFixNode = tempNode;
				// if we are currently merging into the previous branch, then we consider the current branch complete
				// end the current branch search as soon as we find a node to merge into the next branch(therefore, merge into next branch as early as possible)
				// when topMostTendency is 1, this piece of code never excutes, then next round, every 2nd branches will try to merge into its next branch
				// next round, then every 3rd and 2nd branches will try to merge into their corresponding 1st branch. until topMostTendency is entry.size(), where
				// each branch will try to merge into the first branch. This is used to deal with the bottleneck cases and urge branches to merge early
				// the new branch will try to merge into the next branch as early possible, and in the best(shortest) way
				if((entries.size()-entryIndex-1)%topMostTendency!=0)
				{
					if(onNextBranch(tempNode, entryIndex) && !onNextBranch(temp, entryIndex) && earlyMergeNode==NULL)	// this will only be entered once
						earlyMergeNode = tempNode;
				}
			}
		}		// visited all adjacent nodes of the current node
	}
	if(!endingFixNode)		// if we didn't find a fix node, which should NOT happen
		return false;
	else
		setTreeBranchUpForTreeTautening(endingFixNode, start, rnp);			// The branch is complete, trace it back
	return true;
}

// process a tautened tree brach so that each tree node knows which its outgoing edge in the tree is, by marking its treeOutEdgeIndex variable
// the function also helps setting the new rnps of the branch(make sure that the new tree edges have new rnp values)
void RoutingDAG::treeBranchPostProcessingForTreeTautening(Node* start, float rnp)			// the node pointer start always points to an entry node
{
	if(start->getDrawingRNP()==0)								// nothing to process here
		return;
	Node* temp = start;
	while(temp->getNodeType()!=FIX_NODE)
	{
		temp->setDrawingRNP(max(temp->getDrawingRNP(), rnp));	
		// and update the layerUsedIndexReverseDirection vector to be up to date, based on the new branch information
		layerUsedIndexReverseDirection[temp->getLayer()] = min(layerUsedIndexReverseDirection[temp->getLayer()], temp->getLayerIndex());
		for(int i=0; i<temp->getOutSize(); i++)
		{
			if(temp->getOutEdge(i)->treeEdgeOrNot())			// if we find a tree edge
			{
				temp->setTreeOutEdgeIndex(i);
				break;			// break out of the for loop so that we get to the next level in the tree
			}
		}						// then update the drawingRNP information of the edge
		temp->getOutEdge(temp->getTreeOutEdgeIndex())->setDrawingRNP(max(temp->getOutEdge(temp->getTreeOutEdgeIndex())->getDrawingRNP(), rnp));
		temp = temp->getOutNode(temp->getTreeOutEdgeIndex());	// next level of tree node
	}			// after the function, if we have the start node, we could use the treeOutEdgeIndex variable to trace the entire tree branch
	temp->setDrawingRNP(max(temp->getDrawingRNP(), rnp));		// now temp is a fix node, do the same: setting rnp and updating layerUsedIndexReverseDirection
	layerUsedIndexReverseDirection[temp->getLayer()] = min(layerUsedIndexReverseDirection[temp->getLayer()], temp->getLayerIndex());
}

// set all the edges that we used to be in the tree, and it is guaruateed to be a tautend tree branch
// when we find the last node of a tree(a fix node), mark the tree from the current node all way up to the entry node start
void RoutingDAG::setTreeBranchUpForTreeTautening(Node* current, Node* start, float rnp)
{
	while(current!=start)
	{
		// if the current new edge is part of the previous branch, then it means we are NOT adding an indegree on the current node
		if(!(current->treeNodeOrNot() && current->getPrevTreeNode()->treeNodeOrNot()))
			current->addInDegree();										// one more edge is coming in	
		current->setTreeNode();											// the node is a tree node
		if(rnp>current->getDrawingRNP())
			current->setDrawingRNP(rnp);								// the max rnp is updated
		current->getPrevTreeEdge()->setTreeEdge();						// the edge leading to this node is a tree edge
		if(rnp>current->getPrevTreeEdge()->getDrawingRNP())
			current->getPrevTreeEdge()->setDrawingRNP(rnp);
		layerUsedIndex[current->getLayer()] = current->getLayerIndex();	// in current layer, the layerIndex node was the last node used in the tree
		current = current->getPrevTreeNode();
	}
}

// set the drawingRNP values from the merging node to its the fix node, only rnp values are updated because the other status variables are already set 
// to be in the tree. This funtion is only called when we are merging a new branch into an existing branch in the tautened tree
void RoutingDAG::setTreeBranchUpMergingForTreeTautening(Node *start, float rnp)
{
	Node* temp = start;								// find the fix node of the branch containing node start
	while(temp->getNodeType()!=FIX_NODE)
	{
		temp->setDrawingRNP(rnp);
		temp->getOutEdge(temp->getTreeOutEdgeIndex())->setDrawingRNP(rnp);
		temp = temp->getOutNode(temp->getTreeOutEdgeIndex());
	}
	temp->setDrawingRNP(rnp);						// temp is now a fix node
}


/***************************************************************************************************************************************************************/
//HELPER FUNCTIONS USED TO HELP GENERATING TREE BRANCHES

// find the next node of node current on the next branch of the current branch starting from entries[entryIndex], if none, return NULL
Node* RoutingDAG::nextNodeOfGivenNodeOnNextBranch(Node* current, int entryIndex)
{
	int nextBranchEntryIndex = findFeasibleNextEntryNode(entryIndex);
	if(nextBranchEntryIndex==-1)
		return NULL;
	return nextNodeOfGivenNodeOnCurrentBranch(current, nextBranchEntryIndex);
}

// find the next node of node current on the current branch starting from entries[entryIndex], if none, return NULL
Node* RoutingDAG::nextNodeOfGivenNodeOnCurrentBranch(Node *current, int entryIndex)
{
	Node* temp = entries[entryIndex];
	while(temp->getNodeType()!=FIX_NODE)
	{
		if(temp==current)
			return temp->getOutNode(temp->getTreeOutEdgeIndex());
		temp = temp->getOutNode(temp->getTreeOutEdgeIndex());
	}
	return NULL;
}

// given a starting entry node of the current branch, find the merging node of previous 2 branches of the bottommost tree(if there is any)
// example: current starting node in (layer, layerIndex) format is (0, 3). This function finds the merging node of original branches
// starting from (0, 2) and (0, 1). Reason for testing is that we always want to leave enough space for the remaining to be routed branches. If it is 
// a merging node, then we are not going to use it now. Otherwise, the remaining branches will not have a way to get routed.
Node* RoutingDAG::findMergingNodeOfPrevious2BranchesInBottomTree(int entryIndex)
{
	int previousBranch = findFeasiblePreviousEntryNode(entryIndex);
	if(previousBranch==-1)						// if there is NO previous branch existing
		return NULL;
	Node* head1 = entries[previousBranch];		// the starting node of the previous branch
	previousBranch = findFeasiblePreviousEntryNode(previousBranch);
	if(previousBranch==-1)						// if there is NOT a second previous branch
		return NULL;							// then of course current cannot be a merging node
	Node* head2 = entries[previousBranch];		// the starting node of the second previous branch
	// now we find the merging node of branches head1 and head2
	while(head1!=head2)
	{
		if(head1->getLayer()<head2->getLayer())
			head1 = head1->getOutNode(head1->getTreeOutEdgeIndex());
		else head2 = head2->getOutNode(head2->getTreeOutEdgeIndex());
		// if they are both on the fix level and they are still different, then there is NO merging nodes at all 
		if(head1->getLayer()==numLayers-1 && head2->getLayer()==numLayers-1 && head1!=head2)	
			return NULL;
	}
	return head1;								// when head1 and head2 point to the same node, then its the merging node
}

// test the previous branch with the current edge so that till the current layer,there was no collision 
bool RoutingDAG::testPreviousBranchTillCurrentLayerEdge(Edge *current, float rnp, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int temp = findFeasiblePreviousEntryNode(entryIndex);
	if(temp==-1)
		return false;
	// we have the starting node of the adjacent feasible branch now
	return testBranchWithEdge(current, rnp, entries[temp], 2);
}

// test the next branch with the current edge so that till the current layer,there was no collision 
bool RoutingDAG::testnextBranchTillCurrentLayerEdge(Edge* current, float rnp, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int temp = findFeasibleNextEntryNode(entryIndex);
	if(temp==-1)
		return false;
	// we have the starting node of the adjacent feasible branch now
	return testBranchWithEdge(current, rnp, entries[temp], 2);
}

// in the next a couple of functions, the return value is true if the node/edgt collides with the branch, is false if there is no collision
// test the current node with the previous branch, the branch is extracted from layerUsedIndex[0]
bool RoutingDAG::testPreviousBranchNode(Node* current, float rnp, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int temp = findFeasiblePreviousEntryNode(entryIndex);
	if(temp==-1)
		return false;
	// we have the starting node of the adjacent feasible branch now
	return testBranchWithNode(current, rnp, entries[temp]);		
}

// test the current edge with the previous branch, the branch is extracted from layerUsedIndex[0]
bool RoutingDAG::testPreviousBranchEdge(Edge* current, float rnp, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int temp = findFeasiblePreviousEntryNode(entryIndex);
	if(temp==-1)
		return false;
	// we have the starting node of the adjacent feasible branch now
	return testBranchWithEdge(current, rnp, entries[temp], 1);
}

// test the current node with the previous branch, the branch is extracted from layerUsedIndex[0]
bool RoutingDAG::testNextBranchNode(Node* current, float rnp, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int temp = findFeasibleNextEntryNode(entryIndex);
	if(temp==-1)
		return false;
	// we have the starting node of the adjacent feasible branch now
	return testBranchWithNode(current, rnp, entries[temp]);		
}

// test the current edge with the previous branch, the branch is extracted from layerUsedIndex[0]
bool RoutingDAG::testNextBranchEdge(Edge* current, float rnp, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int temp = findFeasibleNextEntryNode(entryIndex);
	if(temp==-1)
		return false;
	// we have the starting node of the adjacent feasible branch now
	return testBranchWithEdge(current, rnp, entries[temp], 1);
}

// test the current node with a branch starting at entry node start(bottom most branch searching, so the new branch cannot intersect with previous one)
// return true if there is an intersection
bool RoutingDAG::testBranchWithNode(Node* current, float rnp, Node* start)
{
	//tracking down the branch starting at the Node start, first test the collision with the startnode
	if(current->collisionWithNode(start, rnp))
		return true;
	// then track down the branch
	Node* temp = start;
	while(temp->getNodeType()!=FIX_NODE)			// traverse till the fix node, but no need to test the fix node (cause it will eventually merge into a fix node)
	{
		if(current->collisionWithEdge(temp->getOutEdge(temp->getTreeOutEdgeIndex()), rnp))
			return true;
		Node* tempNode = temp->getOutNode(temp->getTreeOutEdgeIndex());
		if(current->collisionWithNode(tempNode, rnp))
			return true;
		temp = tempNode;							// going down to the next level
	}
	return false;
}

// test the current edge with a branch starting at entry node start(bottom most branch searching, so the new branch cannot intersect with previous one)
// testType 1: test till the end of the branch. testType 2: test till the current layer of current
bool RoutingDAG::testBranchWithEdge(Edge* current, float rnp, Node* start, int testType)
{
	//tracking down the branch starting at the Node start, first test the collision with the startnode
	if(current->collisionWithNode(start, rnp))
		return true;
	// then track down the branch
	Node* temp = start;
	if(testType==1)
	{
		while(temp->getNodeType()!=FIX_NODE)			// traverse till the fix node, but no need to test the fix node (cause it will eventually merge into a fix node)
		{
			if(current->collisionWithEdge(temp->getOutEdge(temp->getTreeOutEdgeIndex()), rnp))
				return true;
			Node* tempNode = temp->getOutNode(temp->getTreeOutEdgeIndex());
			if(current->collisionWithNode(tempNode, rnp))
				return true;
			temp = tempNode;							// going down to the next level
		}
		return false;
	}
	// when testType is 2, we only test till we reach the current edge's tail's layer
	while(temp->getOutNode(temp->getTreeOutEdgeIndex())->getLayer()<current->getTail()->getLayer())
	{
		if(current->collisionWithEdge(temp->getOutEdge(temp->getTreeOutEdgeIndex()), rnp))
			return true;
		Node* tempNode = temp->getOutNode(temp->getTreeOutEdgeIndex());
		if(current->collisionWithNode(tempNode, rnp))
			return true;
		temp = tempNode;							// going down to the next level
	}
	return false;
}

// helper function to test if a node is on the previous branch of the tree, start is the entry of the current branch we are routing
bool RoutingDAG::onPreviousBranch(Node* current, int entryIndex)
{
	// find the previous entry node that actually has some demand on it
	int tempIndex = findFeasiblePreviousEntryNode(entryIndex);
	if(tempIndex==-1)
		return false;								// there is no previous branch, return false
	// we have the starting node of the adjacent feasible branch now
	return onBranchStartingAt(current, entries[tempIndex]);				// the entry node of the previous branch	
}

// almost the same as the previous function except that we are looking into the next branch
bool RoutingDAG::onNextBranch(Node* current, int entryIndex)
{
	int tempIndex = findFeasibleNextEntryNode(entryIndex);
	if(tempIndex==-1)
		return false;
	return onBranchStartingAt(current, entries[tempIndex]);
}

// find the previous branch that actually has demand input (test if the rnp is 0), return -1 if there isn't any
// current branch starts at node start
int RoutingDAG::findFeasiblePreviousEntryNode(int entryIndex)
{
	if(entryIndex==0)	return -1;						// start is the first entry, so no previous entry node
	// find the previous branch that actually has demand input (test if the rnp is 0)
	int previousFeasibleEntryNodeIndex = entryIndex-1;
	while(previousFeasibleEntryNodeIndex!=-1 && entries[previousFeasibleEntryNodeIndex]->getDrawingRNP()==0)
		previousFeasibleEntryNodeIndex--;
	// we have the starting node of the adjacent feasible branch now
	return previousFeasibleEntryNodeIndex;				// this is -1 if the previous few entry nodes have no demand at all
}

// find the next entry node that actually has demand input, return -1 if there isn't any
int RoutingDAG::findFeasibleNextEntryNode(int entryIndex)
{
	int infeasibleIndex = layerStartingIndex[1];		// the start of the second level(entries are on the first level)
	if(entryIndex==infeasibleIndex-1)	return -1;		// start is the last entry, so no previous entry node
	int nextFeasibleEntryNodeIndex = entryIndex+1;
	while(nextFeasibleEntryNodeIndex!=infeasibleIndex && entries[nextFeasibleEntryNodeIndex]->getDrawingRNP()==0)
		nextFeasibleEntryNodeIndex++;
	if(nextFeasibleEntryNodeIndex==infeasibleIndex)		// this happens if the next few entry nodes have no demand at all
		return -1;
	return nextFeasibleEntryNodeIndex;
}

// test if a node current is on the branch starting at entry node start
bool RoutingDAG::onBranchStartingAt(Node* current, Node* start)
{
	while(start->getNodeType()!=FIX_NODE)
	{
		if(current==start)
			return true;
		start = start->getOutNode(start->getTreeOutEdgeIndex());	// get the next node in the branch
	}
	if(current==start)	return true;				// current is the fix node of the previous branch
	return false;
}

// test to make sure that the minimum distance between 2 merging nodes constraint is satisfied
// return true is the Node current is too close to the branches' merging nodes
bool RoutingDAG::testDistanceTooCloseToMergingNodesOnNextBranch(Node *current, int entryIndex)
{
	if(minimumDistanceBetweenMergingNodes == 0)		// no need to test distance at all
		return false;
	int nextEntryIndex = findFeasibleNextEntryNode(entryIndex);
	if(nextEntryIndex==-1)							// there is no next entry
		return false;
	return testDistanceTooCloseToMergingNodesOnBranch(current, entries[nextEntryIndex]);
}

// test to make sure that the minimum distance between 2 merging nodes constraint is satisfied
// return true is the Node current is too close to the branches' merging nodes
bool RoutingDAG::testDistanceTooCloseToMergingNodesOnPreviousBranch(Node* current, int entryIndex)
{
	if(minimumDistanceBetweenMergingNodes == 0)		// no need to test distance at all
		return false;
	int previousEntryIndex = findFeasiblePreviousEntryNode(entryIndex);
	if(previousEntryIndex==-1)						// there is no previous entry
		return false;
	return testDistanceTooCloseToMergingNodesOnBranch(current, entries[previousEntryIndex]);
}

// test to make sure that the minimum distance between 2 merging nodes constraint is satisfied
// return true is the Node current is too close to the branches' merging nodes
bool RoutingDAG::testDistanceTooCloseToMergingNodesOnBranch(Node *current, Node *start)
{
	while(start)
	{
		if(start!=current && start->getInDegree()==2)	// if we get a merging node here, test the distance
		{											// if the node current is too close to this merging node
			if(sqrt((current->getX()-start->getX())*(current->getX()-start->getX()) +
					(current->getY()-start->getY())*(current->getY()-start->getY()) +
					(current->getZ()-start->getZ())*(current->getZ()-start->getZ()))<minimumDistanceBetweenMergingNodes)
				return true;
		}
		if(start->getNodeType()==FIX_NODE)
			break;
		start = start->getOutNode(start->getTreeOutEdgeIndex());
	}
	return false;									// was never too close, then return NEVER too close
}

// test if a node is too close to the previous merging node on the current branch (if there is one), note that current is on the current branch
bool RoutingDAG::testDistanceTooCloseToMergingNodesOnCurrentBranch(Node* toBeDecided, Node *current, int entryIndex)
{
	// if its impossible to have a merging node on the current branch OR there is NO distance requirement 
	// OR the current node is NOT part of the next branch at all, if so, then there is NO possible merging node along the current branch
	if(minimumDistanceBetweenMergingNodes<=0 || findFeasibleNextEntryNode(entryIndex)==-1 || !onNextBranch(current, entryIndex))
		return false;
	Node* start = entries[entryIndex];			// the start of the current branch. Entry node cannot be a merging node
	while(current!=start)						// walk up the branch to find the nearest merging node if there is one
	{
		Node* prevNodeOnCurrentBranch = current->getPrevTreeNode();	// the previous node on the current branch
		if(!onNextBranch(prevNodeOnCurrentBranch, entryIndex))		// if the previous node is NOT on the next branch, then the current node is a merging node
			break;
		current = prevNodeOnCurrentBranch;							// walk up the current branch
	}
	// test the distance between the merging node and the toBeDecided Node
	if(sqrt((toBeDecided->getX()-current->getX())*(toBeDecided->getX()-current->getX()) + (toBeDecided->getY()-current->getY())*(toBeDecided->getY()-current->getY())
			+(toBeDecided->getZ()-current->getZ())*(toBeDecided->getZ()-current->getZ()))<minimumDistanceBetweenMergingNodes)
		return true;
	return false;
}

/***************************************************************************************************************************************************************/
// after the tree is generated, for each node and edge in the tree, calculate their operational flxibilities
void RoutingDAG::generateOperFlexPairs(float* radii, int length, vector<WeatherData> &wData, float effectiveThres)
{
	if(status==TREE_NOT_GENERATED)
	{
		std::cerr << "\nThe Tree Was Not Generated Yet." << std::endl;
		return;
	}
	// generate the pairs for each branch of the tree
	for(int i=0; i<entries.size(); i++)
	{
		Node* tempNode = entries[i];
		if(entries[i]->getDrawingRNP()==0)			// if there is no branch coming out of the current entry node
			continue;
		while(tempNode->getNodeType()!=FIX_NODE)
		{
			Edge* tempEdge = tempNode->getOutEdge(tempNode->getTreeOutEdgeIndex());
			// we are going to recompute operational flxity pairs, before doing so, clear the vectors
			tempEdge->reset();						// clear the edge's operational flxity related vectors
			tempNode->clearFreeRadiusVector();		// clear the node's operational flxity related vectors		
			// insert operational flexibility pairs of the node and its outgoing tree edge
			for(int j=0; j<length; j++)
			{
				double tempProb = tempNode->testRadiusWithWeatherDataSet(radii[j], wData, effectiveThres);
				tempNode->insertFreeRadiusVec(radii[j], tempProb);
				tempProb = tempEdge->testRNPWithWeatherDataSet(radii[j], wData, effectiveThres);
				tempEdge->insertOperFlex(radii[j], tempProb, 1);		// insert in to the rnp vector
				tempProb = tempEdge->testPathStretchWithWeatherDataSet(radii[j], wData, effectiveThres);
				tempEdge->insertOperFlex(radii[j], tempProb, 2);		// insert in to the path stretching vector
				tempProb = tempEdge->testWiggleRoomWithWeatherDataSet(radii[j], wData, effectiveThres);
				tempEdge->insertOperFlex(radii[j], tempProb, 3);		// insert in to the wiggle room vector
			}
			/***************************************************************************************************************/
			// generate the deviation node for the edge, each edge has at most 3 deviation nodes, 
			// the minimum interval between nodes is preset to be rnp*2; Each node has its own operational flexibility pairs
			double edgeLength = tempEdge->getLength();
			float currentRNP = tempEdge->getDrawingRNP();
			if(edgeLength>=8*currentRNP)										// generate the nodes at 1/4 position from head to tail
			{
				Node* temp = new Node((tempEdge->getHead()->getX()*3+tempEdge->getTail()->getX())/4, (tempEdge->getHead()->getY()*3+tempEdge->getTail()->getY())/4,
									  (tempEdge->getHead()->getZ()*3+tempEdge->getTail()->getZ())/4);
				// compute the operational flxibility pairs of each generated deviation node
				for(int j=0; j<length; j++)
					temp->insertFreeRadiusVec(radii[j], tempNode->testRadiusWithWeatherDataSet(radii[j], wData, effectiveThres));
				tempEdge->insertOperFlexDeviationCandidateNode(temp);
			}
			if(edgeLength>4*currentRNP)										// generate the nodes at 1/2 position from head to tail
			{
				Node* temp = new Node((tempEdge->getHead()->getX()+tempEdge->getTail()->getX())/2, (tempEdge->getHead()->getY()+tempEdge->getTail()->getY())/2,
									  (tempEdge->getHead()->getZ()+tempEdge->getTail()->getZ())/2);
				// compute the operational flxibility pairs of each generated deviation node
				for(int j=0; j<length; j++)
					temp->insertFreeRadiusVec(radii[j], tempNode->testRadiusWithWeatherDataSet(radii[j], wData, effectiveThres));
				tempEdge->insertOperFlexDeviationCandidateNode(temp);
			}
			if(edgeLength>=8*currentRNP)										// generate the nodes at 3/4 position from head to tail
			{
				Node* temp = new Node((tempEdge->getHead()->getX()+3*tempEdge->getTail()->getX())/4, (tempEdge->getHead()->getY()+3*tempEdge->getTail()->getY())/4,
									  (tempEdge->getHead()->getZ()+3*tempEdge->getTail()->getZ())/4);
				// compute the operational flxibility pairs of each generated deviation node
				for(int j=0; j<length; j++)
					temp->insertFreeRadiusVec(radii[j], tempNode->testRadiusWithWeatherDataSet(radii[j], wData, effectiveThres));
				tempEdge->insertOperFlexDeviationCandidateNode(temp);
			}
			/***************************************************************************************************************/
			tempNode = tempNode->getOutNode(tempNode->getTreeOutEdgeIndex());		// walk to the next step of the tree
		}
	}
	// finally, computing the operational flexibility of the fix nodes
	for(int i=0; i<fixes.size(); i++)
	{
		Node* tempNode = fixes[i];
		for(int j=0; j<length; j++)
			tempNode->insertFreeRadiusVec(radii[j], tempNode->testRadiusWithWeatherDataSet(radii[j], wData, effectiveThres));
	}
}
