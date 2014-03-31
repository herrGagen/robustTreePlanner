:: Compile using "MSBuild", then call the necessary commands to generate a tree.

:: Compile
MSBuild.exe /p:Configuration=Release RobustTree.sln

CALL run.bat %*
