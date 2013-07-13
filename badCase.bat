CALL compile_and_run.bat -lanewidth 4 -angle 270 -s 0 -o 0 -dthresh 0.7 -nethresh 0.5 -oname outout1 -twname outout1 -cinput originalInput.txt -iname 0600 -angle 270 -quadrantsize 90 -cellwidth 2.761195 -operflex 1.6567 2.2089 2.7612 3.3134 3.8657 -fixnodes 1 -demandshift 1 -demanddrop 8

RobustTree.exe -cinput Quadrant_4.txt

echo Checking outout1.xml
ruby bugCheck.rb outout1.xml

echo checking Quadrant_4.xml
ruby bugCheck.rb Quadrant_4.xml
