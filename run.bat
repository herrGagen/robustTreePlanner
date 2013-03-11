:: This file is identical to "compile_and_run.bat" except it skips the compile step.

@echo OFF
@echo %time%

:: Convert weather to a format useable by RTP, pass it all the command line parameters
:: This also calls ``create_input.rb`` to create the "inputs.txt" file
ruby convert_files.rb %*

:: Run the main program
.\RobustTree.exe

@echo %time%