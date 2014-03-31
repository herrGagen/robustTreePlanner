#ifndef DYNPROGTREEGENERATORTESTER_H
#define DYNPROGTREEGENERATORTESTER_H

#include "DynProgTreeGenerator.h"
#include "UserInterface.h"
#include "gtest/gtest.h"

class DynProgTreeGeneratorTester
{
public:
	class DynProgTreeGenerator generator;
public:
	DynProgTreeGeneratorTester( unsigned int NUMR = 6, unsigned int NUMDEM = 6);
	DynProgTreeGeneratorTester(  std::string inputFileName );

	bool indexFailsOnZeroDemands() const;
	bool indexCountsSequentially() const;
	bool reverseIndexVerifier() const;
	unsigned int numBranchesNotTerminatingAtArrivalNodes();
	std::pair<double,double> getInputAndOutput() ;
};

class DynProgTreeGeneratorTester_Default : public ::testing::Test, public DynProgTreeGeneratorTester
{
public:
	DynProgTreeGeneratorTester_Default(): DynProgTreeGeneratorTester(6,6){}
};

class DynProgTreeGeneratorTester_TXT : public ::testing::Test, public DynProgTreeGeneratorTester
{
public:
	DynProgTreeGeneratorTester_TXT(): DynProgTreeGeneratorTester("six.txt"){}
};


#endif // already included