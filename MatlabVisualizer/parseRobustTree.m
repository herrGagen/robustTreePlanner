
function [nodes,arcs,branch,weather] = parseRobustTree(a)

nodes = parseNodes(a.Children(1));
arcs = parseArcs(a.Children(2));
branch = parseBranches(a.Children(3));
if(length(a.Children) > 3 )    
    weather = parseXMLWeather(a.Children(4) );
else
    weather = [];
end
