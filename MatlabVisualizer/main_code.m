function visualizeXML( inputFileName )
clear; clc; close all

plot_nas(1)
% parse robust tree
if( nargin == 0)
    [FileName,PathName,FilterIndex] = uigetfile('./*.xml')
end
route = readRobustTree(inputFileName)

% load parsed weather file
load wx

% select wx
i = 3;
% 1 - 060000
% 2 - 061500
% 3 - 063000
% 4 - 064500
% 5 - 070000
% 6 - 071500
% 7 - 073000
% 8 - 074500
% 9 - 080000

% specify altitude slice
alt_low = 8000;
alt_high = 26000;
% specify dev prob
dev_thld = 0.75;

% select ensemble member to display
%n = [1 2 3 4 5 6];
n = [1];
% plot ensemble members
for j = n
    lat = wx(i).member(j).dat(:,1);
    lon = wx(i).member(j).dat(:,2);
    alt = wx(i).member(j).dat(:,3);
    dev = wx(i).member(j).dat(:,4);    
    k = find(alt>=alt_low & alt<=alt_high & dev>=dev_thld);    
    plot(lon(k),lat(k),'r.')    
end