
function varargout = dirs(varargin);

if ~isempty(varargin)
    f = dir(varargin{1});
else
    f = dir;
end
f = f(3:end);

for i = 1:length(f)
    d(i) = datenum(f(i).date);
    f(i).date = datestr(d(i),'mm-dd-yyyy  HH:MM AM');
end
[~,b] = sort(d);
f = f(b);

if nargout == 1
    varargout{1} = f;
else
    varargout = [];
    disp(f)
end

end

%%
function disp(f)
for i = 1:length(f)
    fprintf('%s \t %s\n',f(i).date,f(i).name)
end
end