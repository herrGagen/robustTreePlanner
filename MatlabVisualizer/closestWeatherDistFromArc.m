function minDist = closestWeatherDistFromArc(arcLat1, arcLon1, arcLat2, arcLon2, weatherLat, weatherLon)
%
% The distance from a line segment (x1,y1)-(x2,x2) to a point (a,b) is
% found by first parameterizing the line segment to:
% (x,y) = (x1,y1) + t(x2-x1,y2-y1) = (x1,y1) + t(dx,dy)
% Then we take the distance to point (a,b):
% d^2 = ( a-x1 - tdx )^2 + ( b-y1 - tdy )^2
% This is minimized when the partial derivative wrt t is 0
% 0 = -2dx( a-x1 - tdx ) -2dy( b - y1 - tdy )
% 2dx(a-x1) + 2dy(b-y1) = 2tdx^2 + 2tdy^2
% t = ( dx(a-x1) + dy(b-y1) )/ (dx^2 + dy^2)

% Because we are working in lat/lon rather than x,y space, we
% want to convert lat lon to distances. We're going to locally approximate
% Distance as linear in lat/lon using the following 2 ocnversions

appxNmPerLat = haversine([arcLat1 arcLon1],[arcLat2 arcLon1])/abs(arcLat1-arcLat2);
appxNmPerLon = haversine([arcLat1 arcLon1],[arcLat1 arcLon2])/abs(arcLon1-arcLon2);

% This changes our distance calculation to:
% d^2 = appxNmPerLat^2( a-x1 - tdx )^2 + appxNmPerLon^2( b-y1 - tdy )^2
% And our resulting minimum to:
% t = ( appxNmPerLat^2*dx(a-x1) + appxNmPerLon^2*dy(b-y1) ) /
%     (appxNmPerLat^2*dx^2 + appxNmPerLon^2*dy^2)

dLat = arcLat2 - arcLat1;
dLon = arcLon2 - arcLon1;

dists(1,:) = distanceFromPointOnArcToWeather(0, arcLat1, arcLon1, dLat, dLon, weatherLat(:), weatherLon(:), appxNmPerLat, appxNmPerLon);
dists(2,:) = distanceFromPointOnArcToWeather(1, arcLat1, arcLon1, dLat, dLon, weatherLat(:), weatherLon(:), appxNmPerLat, appxNmPerLon);

num = appxNmPerLat^2*dLat*(weatherLat-arcLat1) +  ...
      appxNmPerLon^2*dLon*(weatherLon-arcLon1) ;
denom = ( appxNmPerLat^2*dLat^2 + appxNmPerLon^2*dLat^2 );

closestT = num/denom;

dists(3,:) = dists(2,:);
goodInds = (closestT(:) > 0) & (closestT(:) < 1);
dists(3,goodInds) = distanceFromPointOnArcToWeather(1, arcLat1, arcLon1, dLat, dLon, weatherLat(goodInds), weatherLon(goodInds), appxNmPerLat, appxNmPerLon);

minDist = min(dists(:) );

end %function distanceFromArcToPoint

function dist = distanceFromPointOnArcToWeather(t, arcLat1, arcLon1, dLat, dLon, weatherLat, weatherLon, appxNmPerLat, appxNmPerLon)
% dist = distanceFromPointOnArcToPoint(t, arcLat1, arcLon1, dLat, dLon, weatherLat, weatherLon, appxNmPerLat, appxNmPerLon)
%
% Returns the distance from arc (arcLat1,arcLat2) + t(dLat,dLon) to point
% (weatherLat, weatherLon) 


dist = zeros(size(weatherLat) );
for i = 1:length(weatherLat)
    dist(i) = haversine([arcLat1 + t*dLat arcLon1 + t*dLon], [weatherLat(i) weatherLon(i)]);
end

end