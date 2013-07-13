
function branch = parseBranches(nodes)

branch = [];
for i = 1:length(nodes.Children)
    a = nodes.Children(i);
    b = [];
    for j = 1:length(a.Children)        
        b(j) = str2num(a.Children(j).Attributes.Value);
    end
    branch(i).ind = b;
end