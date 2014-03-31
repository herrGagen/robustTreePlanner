
function [crdC] = getCenterLine(arcs,nodes)

ind = [nodes.index];
latC = []; lonC = []; indC = [];
for i = 1:length(arcs)
    k0 = find(ind==arcs(i).start_node);
    latC = [latC nodes(k0).latitude];
    lonC = [lonC nodes(k0).longitude];
    indC = [indC arcs(i).start_node];
end
k1 = find(ind==arcs(i).end_node);
latC = [latC nodes(k1).latitude];
lonC = [lonC nodes(k1).longitude];
indC = [indC arcs(i).end_node];

crdC = [latC; lonC; indC]';