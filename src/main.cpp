#include "UserInterface.h"
#include <iostream>

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
  
	UserInterface* userInterface = new UserInterface();
	userInterface->ProgramBegins(inputFile);
	delete userInterface;
	return 0;
}