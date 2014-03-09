function minDist = closestWeatherDistFromArc(arcLat1, arcLon1, arcLat2, arcLon2, weatherLat, weatherLon)
% minDist = closestWeatherDistFromArc(arcLat1, arcLon1, arcLat2, arcLon2, weatherLat, weatherLon)
%
% Using haversine distance, finds the distance to the closest point in the 
% weather ensemble
% 
% NOTE: If the arc lengths are smaller than the distane to the closest point,
% then this routine will simply report arc lengths.
%
% SEE ALSO: haversine

arcLength = haversine([arcLat1 arcLon1],[arcLat2 arcLon2]);

dLat = arcLat2 - arcLat1;
dLon = arcLon2 - arcLon1;

distLeft = distanceFromPointOnArcToWeather(0, arcLat1, arcLon1, dLat, dLon, weatherLat(:), weatherLon(:));
distRight = distanceFromPointOnArcToWeather(1, arcLat1, arcLon1, dLat, dLon, weatherLat(:), weatherLon(:));

dists = arcLength*ones(size(weatherLat) );
for i = 1:length(weatherLat)
    if(distLeft(i) < arcLength || distRight(i) < arcLength)
        [~, dists(i)] = fminbnd(@(t) distanceFromPointOnArcToWeather(t, arcLat1, arcLon1, dLat, dLon, weatherLat(i), weatherLon(i) ), 0, 1);
        %[~, dists(i)] = findGlobalMin(@(t) distanceFromPointOnArcToWeather(t, arcLat1, arcLon1, dLat, dLon, weatherLat(i), weatherLon(i) ), 0, 1);
    end
end

minDist = min( dists(:) );

end %function distanceFromArcToPoint

function dist = distanceFromPointOnArcToWeather(t, arcLat1, arcLon1, dLat, dLon, weatherLat, weatherLon)
% dist = distanceFromPointOnArcToPoint(t, arcLat1, arcLon1, dLat, dLon, weatherLat, weatherLon)
%
% Returns the distance from arc (arcLat1,arcLat2) + t(dLat,dLon) to point
% (weatherLat, weatherLon) 

if(length(t) ~= length(weatherLat) )
    t = t*ones(size(weatherLat) );
end
dist = zeros(size(weatherLat) );
for i = 1:length(weatherLat)
    dist(i) = haversine([arcLat1 + t(i)*dLat arcLon1 + t(i)*dLon], [weatherLat(i) weatherLon(i)]);
end

end