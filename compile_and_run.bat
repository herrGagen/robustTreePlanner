:: Compile using "MSBuild", then call the necessary commands to generate a tree.

:: Compile
MSBuild.exe RobustTree.sln

:: Convert weather to a format useable by RTP, pass it all the command line parameters
:: This also calls ``create_input.rb`` to create the "inputs.txt" file
ruby convert_files.rb %*

:: Run the main program
.\RobustTree.exe %*