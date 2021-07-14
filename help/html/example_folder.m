function example_folder()
%find the location of the asmita demo folder and open it
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'CoastalTools'));
fpath = [appinfo(idx(1)).location,'/example'];
open(fpath)