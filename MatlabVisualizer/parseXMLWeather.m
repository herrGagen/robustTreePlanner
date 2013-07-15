function weather = parseXMLWeather(nodes)

latVec = zeros(size(nodes.Children));
lonVec = zeros(size(nodes.Children));
for i = 1:length(nodes.Children)
   latVec(i) = str2double( nodes.Children(i).Attributes(1).Value );
   lonVec(i) = str2double( nodes.Children(i).Attributes(2).Value );
end

weather.lat = latVec;
weather.lon = lonVec;