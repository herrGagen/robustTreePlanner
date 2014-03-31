function wx = parseWeather(weatherDir)

wx = [];
startingDir = pwd;
cd(weatherDir);
fs = dirs;
for i = 1:length(fs)
    f = fs(i).name    
    if ~isempty(strfind(f,'.txt'))
        [~,~,~,lat,lon,alt,dev] = textread(f,'%f %f %f %f %f %f %f','delimiter',',','headerlines',11);
        str0 = regexp(f,'[_.]+','split');
        str1 = str0{4};
        n = str2num(str1(7:end));
        str2 = regexp(str0{3},'T','split');
        flg = 0;
        for k = 1:length(wx)
            if strcmp(str2{2},wx(k).id)
                flg = 1;
                break;
            end
        end
        if ~flg
            k = length(wx)+1;
            wx(k).id = str2{2};
        end
        wx(k).member(n).dat = [lat lon alt dev];
    end
end

% save to .mat file
cd(startingDir);