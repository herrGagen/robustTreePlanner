#include <iostream>

#include "UserInterface.h"
#include "TreeVerifier.h"
#include "DynProgTreeGenerator.h"

int main(int argc, char* argv[])
{	

  std::string inputFile = "inputs.txt";
  if (argc > 1)
  {
    for (int i = 0; i < argc; i++)
    {
      std::string argvi = argv[i];
      std::cout << argvi << std::endl;
      if (argvi == "-cinput")
      {
        inputFile = argv[i+1];
      }
    }
  }

  UserInterface userInterface;
	if( userInterface.ProgramBegins(inputFile) == false )
  {
    std::ofstream emptyFile;
    emptyFile.open( userInterface.getOutputFileName() + ".xml");
    emptyFile << "" << std::endl;
    emptyFile.close();
    emptyFile.open( userInterface.getOutputFileName() + ".dp.xml");
    emptyFile << "" << std::endl;
    emptyFile.close();
    return -1;
  }
	TreeVerifier checker(userInterface);
	checker.appendReportToFiles("InvalidTrees.txt","ValidTrees.txt");
  checker.outputTreeLengthStats();
	DynProgTreeGenerator dp(userInterface);
	dp.writeBestTreeToDAG(userInterface);
	TreeVerifier dpChecker(userInterface);
	dpChecker.appendReportToFiles("dpInvalidTrees.txt","dpValidTrees.txt");
    dpChecker.outputTreeLengthStats();
	userInterface.setOutputFileName( userInterface.getOutputFileName() + ".dp.xml" );
	userInterface.saveTreeInformation();
	return 0;
}
