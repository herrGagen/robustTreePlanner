
function nodes = parseNodes(Nodes)
nodes = [];
for i = 1:length(Nodes.Children)
    a = Nodes.Children(i);
    nd = [];
    for j = 1:length(a.Attributes)
        x = str2num(a.Attributes(j).Value);
        if ~isempty(x)
            nd = setfield(nd,a.Attributes(j).Name,x);
        else
            nd = setfield(nd,a.Attributes(j).Name,a.Attributes(j).Value);
        end
    end
    disc = [];
    b = a.Children;
    if ~isempty(b.Children)
        for j = 1:length(b.Children)
            d = [];
            for k = 1:length(b.Children(j).Attributes)
                x = str2num(b.Children(j).Attributes(k).Value);
                if ~isempty(x)
                    d = setfield(d,b.Children(j).Attributes(k).Name,x);
                else
                    d = setfield(d,b.Children(j).Attributes(k).Name,b.Children(j).Attributes(k).Value);
                end
            end
            disc = [disc d];
        end
    end
    nd.disc = disc;
    nodes = [nodes nd];
end
end