#include "UserInterface.h"

int main()
{
	UserInterface* userInterface = new UserInterface();
	userInterface->ProgramBegins();
	delete userInterface;
	return 0;
}