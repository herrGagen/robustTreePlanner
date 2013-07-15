del outout1.xml
del Quadrant_4.xml
del Synthetic-Angle000.xml
del Synthetic-Angle090.xml
del Synthetic-Angle180.xml
del Synthetic-Angle270.xml

CALL compile_and_run.bat -lanewidth 4 -angle 270 -s 0 -o 0 -dthresh 0.7 -nethresh 0.5 -oname outout1 -twname outout1 -cinput originalInput.txt -iname 0600 -angle 270 -quadrantsize 90 -cellwidth 2.761195 -operflex 1.6567 2.2089 2.7612 3.3134 3.8657 -fixnodes 1 -demandshift 1 -demanddrop 8

RobustTree.exe -cinput Quadrant_4.txt

RobustTree.exe -cinput Synthetic-Angle000.txt
RobustTree.exe -cinput Synthetic-Angle090.txt
RobustTree.exe -cinput Synthetic-Angle180.txt
RobustTree.exe -cinput Synthetic-Angle270.txt


