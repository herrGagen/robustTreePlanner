
function plot_nas(varargin)

if ~isempty(varargin)
    isLine = varargin{1};
else
    isLine = 0;
end
no_tcn = 0;
if length(varargin)==2
    no_tcn = 1;
end

% figure
hold on
% set(gca,'FontSize',13)

% set flag to display centers and tracons
centers = 1;
tracons = 1;

% set RGB values
USfill = [0.7529 0.7529 0.7529];
USedge = [0.6863 0.6863 0.6863];
CTedge = [0 0.6000 0.8000];

% limits for displaying nas
xmin = -130; xmax = -65;
ymin = 20; ymax = 55;

% degree to radian
d2r=pi/180;

% determine whether to use map projection
proj = 0;
if proj
    daspect([1 1 1])
else
    daspect([1 cos(37*d2r) 1])
end

%Define names of Tracons for plotting
tcns = {'PCT';'N90';'C90';'SCT';'PHL';'CVG';'ATL'};
% tcns = {'C90'};
% plot US
fid = fopen('usmap_48.dat','r');
file = fscanf(fid,'%f %f',[2 inf]);
file = file';
fclose(fid);

indx = find(file(:,1) == 0);
[row,col] = size(indx);
for i=1:row-1
    lat = file(indx(i)+1:indx(i+1)-1,1);
    lon = file(indx(i)+1:indx(i+1)-1,2);
    if proj
        [x y] = lambc([lat lon]*d2r);
    else
        x = lon;
        y = lat;
    end
    if ~isLine
        fill(x, y, USfill, 'EdgeColor', USedge,'tag','usa')
        %         plot(x, y, '-', 'color',[1 0 1],'linewidth',3,'tag','usaedge','visible','on')
    else
        plot(x, y, '-', 'color',USedge)
    end
end

% Plot Centers
str = textread('centers1.dat','%s');
c = char(str);

%returns rows containing center names
index = find(c(:,1)=='.');

%add one more index
index = [index; length(c)+1];

%go through each center
for i = 1:length(index)-1
    tag = 'ctr';
    if centers == 0
        state = 'off';
    else
        state = 'on';
    end
    C = CTedge;
    ctr = c(index(i),2:4);
    t = find(strcmp(ctr,tcns) == 1);
    if isempty(t) == 0
        %This center is a tracon cluster
        C=[0.4, 0.4, 0.4];
        C=[0.1 0.1 0.1];
        if tracons == 0
            state='off';
        else
            state='on';
        end
        tag='tcn';
    end
    
    num = str2num(c(index(i)+3:index(i+1)-1,:));
    n = length(num);
    m = 2*(1:n/2);
    lat = num(m-1);
    lon = -num(m);
    %ensure first and final points are equal
    if or(lat(1)~=lat(end),lon(1)~=lon(end))
        lat=[lat; lat(1)];
        lon=[lon; lon(1)];
    end
    %projection
    if proj
        [x y] = lambc([lat lon]*d2r);
    else
        x = lon;
        y = lat;
    end
    if isempty(t) || (~isempty(t) && no_tcn==0)
        h = plot(x, y, 'm-','linewidth', 1, 'color', C,'tag',tag);
        set(h,'visible',state,'clipping','on');
        %     if isempty(strfind(c(index(i),2:7),'Z'))
        h = text(mean(x),mean(y),c(index(i),2:7),...
            'Color',C,...
            'horizontalalignment','center',...
            'tag',tag,...
            'fontsize',8);
        %     end
        set(h,'visible',state,'clipping','on');
    end    
end

xlabel('Longitude')
ylabel('Latitude')
axis([xmin xmax ymin ymax])