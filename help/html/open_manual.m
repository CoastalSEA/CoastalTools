function open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'CoastalTools'));
fpath = [appinfo(idx(1)).location,'/doc/CoastalTools manual.pdf'];
open(fpath)
