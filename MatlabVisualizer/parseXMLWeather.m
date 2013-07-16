function weather = parseXMLWeather(nodes)

latVec = zeros(size(nodes.Children));
lonVec = zeros(size(nodes.Children));
idVec = zeros(size(nodes.Children));
probVec = zeros(size(nodes.Children));
for i = 1:length(nodes.Children)
   idVec(i)    = str2double( nodes.Children(i).Attributes(1).Value );
   probVec(i)  = str2double( nodes.Children(i).Attributes(2).Value );
   latVec(i)   = str2double( nodes.Children(i).Attributes(3).Value );
   lonVec(i)   = str2double( nodes.Children(i).Attributes(4).Value );
end

weather.lat = latVec;
weather.lon = lonVec;
weather.ensembleId = idVec;
weather.ensembleProb = probVec;