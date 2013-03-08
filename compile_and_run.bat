
:: Compile
MSBuild.exe RobustTree.sln

:: Convert weather to a format useable by RTP
:: This also calls ``create_input.rb`` to create the "inputs.txt" file
ruby convert_files.rb

:: Run the main program
.\RobustTree.exe