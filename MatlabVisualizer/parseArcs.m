
function arcs = parseArcs(Nodes)

arcs = [];
for i = 1:length(Nodes.Children)
    a = Nodes.Children(i);
    arc = [];
    for j = 1:length(a.Attributes)
        x = str2num(a.Attributes(j).Value);
        if ~isempty(x)
            arc = setfield(arc,a.Attributes(j).Name,x);
        else
            arc = setfield(arc,a.Attributes(j).Name,a.Attributes(j).Value);
        end
    end
    rect = [];
    for j = 1:length(a.Children)
        if strfind(a.Children(j).Name,'RNP')==1
            rnps = a.Children(j);
            for k = 1:length(rnps.Children)
                rnp_prob(k,:) = [str2num(rnps.Children(k).Attributes(1).Value) ...
                    str2num(rnps.Children(k).Attributes(2).Value)];
            end
        elseif strfind(a.Children(j).Name,'Operational')==1
            recs = a.Children(j);
            ofa = [];
            n = 0;
            for k = 1:2:length(recs.Children)
                n = n+1;
                att = recs.Children(k);
                for m = 1:length(att.Attributes)
                    rect = setfield(rect,att.Attributes(m).Name,att.Attributes(m).Value);
                end
                rec = recs.Children(k+1);
                for m = 1:length(rec.Children)
                    crd(m,:) = [str2num(rec.Children(m).Attributes(2).Value) ...
                        str2num(rec.Children(m).Attributes(3).Value)];
                end
                rect = setfield(rect,'coord',crd);
                ofa(n).rect = rect;
            end
            
        end
    end
    arc = setfield(arc,'rnp_prob',rnp_prob);
    arc = setfield(arc,'ofa',ofa);
    arcs = [arcs arc];
end
end