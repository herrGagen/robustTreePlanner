function route = readRobustTree(fname)

fprintf('Reading robust tree from %s',fname)
linearizeXML(fname);

a = parseXML(fname);

[nodes,arcs,branch,weather] = parseRobustTree(a);

k = [arcs.index];

for m = 1:length(branch)
    ind = branch(m).ind;
    k1 = [];
    for i = 1:length(ind)
        k1(i) = find(k==ind(i));
    end
    arcs1 = arcs(k1);
    crdC = getCenterLine(arcs1,nodes);
    route(m).lat = crdC(:,1);
    route(m).lon = crdC(:,2);
    for i = 1:length(arcs1)
        for j = 1:length(arcs1(i).ofa)
            lt = arcs1(i).ofa(j).rect.coord(:,1);
            ln = arcs1(i).ofa(j).rect.coord(:,2);
            seg(i).ofa(j).lat = lt;
            seg(i).ofa(j).lon = ln;
            seg(i).ofa(j).width = arcs1(i).ofa(j).rect.width;
            seg(i).ofa(j).prob = arcs1(i).ofa(j).rect.probability;
            plot([ln; ln(1)],[lt; lt(1)],'-g')
        end
    end
    route(m).seg = seg;    
    plot(crdC(:,2),crdC(:,1),'-b.')
    text(crdC(1,2),crdC(1,1),num2str(m))
end

if(~isempty(weather) )
    plot(weather.lon,weather.lat,'r.');
end
fprintf('... Done\n')