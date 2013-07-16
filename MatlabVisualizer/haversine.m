function nmi = haversine(loc1, loc2)
% HAVERSINE     Compute distance between locations using Haversine formula
%   NMI = HAVERSINE(LOC1, LOC2) returns the distance in nautical miles 
%   between locations LOC1 and LOC2 using the Haversine formula.  LOC1 and 
%   LOC2 are latitude and longitude coordinates expressed as numeric 
%   arrays of decimal degrees (where negative indicates West/South).
%
%
%   Examples
%       haversine([53.1472 -1.8494], [52.2044 0.1406]) returns 170.2563
%
%   Inputs
%       LOC must be a 2-valued numeric array specifying the
%       location in decimal degrees.  
%
%   Notes
%       The Haversine formula is used to calculate the great-circle
%       distance between two points, which is the shortest distance over
%       the earth's surface.
%
%       This program was created using equations found on the website
%       http://www.movable-type.co.uk/scripts/latlong.html

% Created by Josiah Renfree
% May 27, 2010

% Convert all decimal degrees to radians
loc1 = deg2rad(loc1);
loc2 = deg2rad(loc2);

%% Begin calculation

R = 3443.89849;                       % Earth's radius in nmi
delta_lat = loc2(1) - loc1(1);        % difference in latitude
delta_lon = loc2(2) - loc1(2);        % difference in longitude
a = sin(delta_lat/2)^2 + cos(loc1(1)) * cos(loc2(1)) * ...
    sin(delta_lon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
nmi = R * c;                                 % distance in nmi