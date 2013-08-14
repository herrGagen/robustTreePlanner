function visualizeXML( inputFileName )

h = plot_nas(1);
% parse robust tree

if( nargin == 0)
    [inputFileName path] = uigetfile('./XML/*.xml');
    route = readRobustTree(fullfile([path inputFileName]) );
else
    route = readRobustTree(inputFileName);
end

[~, name, ~] = fileparts(inputFileName);
figure(h);
title(name);
saveas(h,name,'fig' );