
function linearizeXML(fname)
fid = fopen(fname);
n = 0;
while ~feof(fid)
    n = n+1;
    str{n} = fgetl(fid);
end
fclose(fid);

fid = fopen(fname,'w');
for i = 1:length(str)
    s = regexprep(str{i},'[\r\t]+','');
    fprintf(fid,'%s',s);
end
fclose(fid);
end