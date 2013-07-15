function minDist = verifyRobustTree(fname)

fprintf('Minimum weather clearance in %s is ',fname)
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
end

minDist = Inf;
for i = 1:length(route)
    latVec = route(i).lat;
    lonVec = route(i).lon;
    for j = 1:(length(latVec)-1)
        thisMinDist = closestWeatherDistFromArc(latVec(j),lonVec(j),latVec(j+1),lonVec(j+1),weather.lat,weather.lon);
        if(thisMinDist < minDist)
            minDist = thisMinDist;
        end
    end
end

fprintf('%f\n',minDist);
outf = fopen('xmlClearances.txt','a');
fprintf(outf,'Minimum weather clearance in %s is %f\n',fname,minDist);
fclose(outf);