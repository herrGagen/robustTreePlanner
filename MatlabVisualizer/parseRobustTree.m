
function [nodes,arcs,branch] = parseRobustTree(a)

nodes = parseNodes(a.Children(1));
arcs = parseArcs(a.Children(2));
branch = parseBranches(a.Children(3));