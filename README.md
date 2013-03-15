LinuxRobustTree
===============

Robust Tree Planner for Linux

To run the tree builder, you need to have a folder called "CWAMEnsembles" in the main root directory (where the .sln and .exe files are).
This should be in the standard format and are the weather files you want analyzed.

Next, you can call either:
compile_and_run.bat or run.bat (they do exactly what their name implies).

I'm not sure when you need to use compile_and_run, and it has the dependency of MSBuild.exe, which should be in your path. (for me that means I must use Windows Powershell to call compile_and_run).

You pass either of these .bat files the command line parameters, so here are the parameters and their meanings:

-s
You call -s to indicate which weather file time you want to start with.  0 indicates the lowest, and since there are normally 90 minutes or so of weather data, you can pick it to be 0, 15, 30, 45, 60, 75, etc.
To call it simply type "-s 15" or your desired starting time (offset from the lowest time).

-o
This is the offset, or time gap that you want to include.
So if you have -s 0 -o 15, you will include the first two weather times.
If you have -s 1 -o 15, you'll include the second.
-s 15 -o 0 would do the same.

Both -o and -s default to zero.

-dthresh
This is the deviation threshold parameter for Shang's tree builder. A higher number means fewer edges are considered, a lower number means many edges are considered.

-nethresh
This is the node-edge threshold. If an edge passes through really bad weather, it will not be considered part of a viable tree.
Again, higher is more restrictive.

These default to 0.8 (and you should write the 0 in front of a decimal, so 0.7, not .7)

For file management, there are the following:
-oname
This is just the name for your output file. It defaults to a timestamp from some point while the program is running if you don't specify anything.
The default will look like this:
2013-03-08_01-49-22.xml, indicating year-month-day_hour-minutes-seconds.xml

-iname
WARNING: IF YOU USE THIS PARAMETER, THE DIRECTORY MUST EXIST AND CONTAIN WEATHER DATA OR THE CODE WILL HAVE AN ERROR
This parameter allows you to select what directory contains the weather data.
The default value would be equivalent to calling
> \run.bat -iname CWAMEnsembles

-twname
This allows you to rename the directory where the temporary weather is stored. By default it gives a string indicating range of times (of weather) it includes.

-cinput
This allows you to specify the name of the file used for the C++ code for its input. The value defaults to "inputs.txt"
(for future maintainers, that value, "inputs.txt" is hardcoded in two places: once in create_input.rb:9, and once in src/main.cpp:5)
