function numViolatedEdges = verifyRobustTree(fname, thresholdDist, edgeProbThresh)
% numViolatedEdges = verifyRobustTree(fname, thresholdDist, edgeProbThresh)

fprintf('Number of bad edges in %s is ',fname)
linearizeXML(fname);

a = parseXML(fname);

[nodes,arcs,branch,weather] = parseRobustTree(a);

k = [arcs.index];

routes( length(branch) ).lat = [];
routes( length(branch) ).lon = [];
for m = 1:length(branch)
    ind = branch(m).ind;
    k1 = [];
    for i = 1:length(ind)
        k1(i) = find(k==ind(i));
    end
    arcs1 = arcs(k1);
    crdC = getCenterLine(arcs1,nodes);
    routes(m).lat = crdC(:,1);
    routes(m).lon = crdC(:,2);
end

% Count hte number of edges that are violated
% Where an edge is defined as being violated if it is closer than
% thresholdDist from a number of weather cells whose emsembleProbabilities
% sum to edgeProbThresh or greater
ensembleIndices = unique(weather.ensembleId);
numViolatedEdges = 0;
for i = 1:length(routes)
    latVec = routes(i).lat;
    lonVec = routes(i).lon;
    for j = 1:(length(latVec)-1)      
        edgeProb = 1;
        for ensemble = ensembleIndices
            ndxInThisEnsemble = (weather.ensembleId == ensemble);
            thisMinDist = closestWeatherDistFromArc(latVec(j),lonVec(j),latVec(j+1),lonVec(j+1),weather.lat(ndxInThisEnsemble),weather.lon(ndxInThisEnsemble) );
            if(thisMinDist < thresholdDist)
                edgeProb = edgeProb - max(weather.ensembleProb(ndxInThisEnsemble) );
            end
        end
        if( edgeProb < edgeProbThresh )
           numViolatedEdges = numViolatedEdges+1;
        end
    end
end

fprintf('%d\n',numViolatedEdges);
outf = fopen('xmlViolations.txt','a');
fprintf(outf,'Number of bad edges in %s is %d\n',fname,numViolatedEdges);
fclose(outf);