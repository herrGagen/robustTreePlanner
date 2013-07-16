function numViolatedFiles = verifyAllXML
% Runs all XML files in the current directory through verifyRobustTree
% and reports how many of them had illegal edges.
%
% This has a false positive rate of 0, but misses cases where the weather
% is between [laneWidth laneWidth+cellSize]

tic;
h = waitbar(0,'Please wait');
numViolatedFiles = 0;

lwNeg = dir('*laneWidth-1*.xml');
lwTwo = dir('*laneWidth2*.xml');
lwFour = dir('*laneWidth4*.xml');
allFiles = dir('*.xml');

numComplete = 0;
numToComplete = length(lwNeg) + length(lwTwo) + length(lwFour) + length(allFiles);

for i = 1:length(lwNeg)
    fprintf('(1/3) %0.5d ',i); 
    temp = verifyRobustTree( lwNeg(i).name, .82, .8); 
    if(temp > 0)
        numViolatedFiles = numViolatedFiles+1;
    end
    numComplete = numComplete + 1;
    waitbar(numComplete/numToComplete,h,timeRemainingMessage(numComplete,numToComplete) );
end

for i = 1:length(lwTwo)
    fprintf('(2/3) %0.5d ',i); 
    temp = verifyRobustTree( lwTwo(i).name, 2, .8); 
    if(temp > 0)
        numViolatedFiles = numViolatedFiles+1;
    end
    numComplete = numComplete + 1;
    waitbar(numComplete/numToComplete,h,timeRemainingMessage(numComplete,numToComplete) );
end

for i = 1:length(lwFour)
    fprintf('(3/3) %0.5d ',i); 
    temp = verifyRobustTree( lwFour(i).name, 4, .8); 
    if(temp > 0)
        numViolatedFiles = numViolatedFiles+1;
    end
    numComplete = numComplete + 1;
    waitbar(numComplete/numToComplete,h,timeRemainingMessage(numComplete,numToComplete) );
end

for i = 1:length(allFiles)
    fprintf('%0.5d ',i); 
    visualizeXML( allFiles(i).name); 
    numComplete = numComplete + 1;
    waitbar(numComplete/numToComplete,h,timeRemainingMessage(numComplete,numToComplete) );
end

for i = 1:length(allFiles)
    fprintf('%0.5d ',i); 
    visualizeXML( allFiles(i).name); 
    numComplete = numComplete + 1;
    waitbar(numComplete/numToComplete,h,timeRemainingMessage(numComplete,numToComplete) );
end


function msg = timeRemainingMessage(numComplete,numToComplete)

duration = toc;
if(numComplete > 0 || duration == 0)
    durPerCall = duration/numComplete;
    totalTime = numToComplete*durPerCall;
    ttString = datestr(totalTime/(24*3600),'HH:MM:SS');
    remainingTime = totalTime - duration;
    rtString = datestr(remainingTime/(24*3600),'HH:MM:SS');
    msg = sprintf('%s remaining of %s',rtString,ttString);
else
    msg = 'No idea how long this will take. Pack a lunch.';
end

