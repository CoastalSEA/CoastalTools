function ct_data_cleanup(muicat,src)
%
%-------function help------------------------------------------------------
% NAME
%   ct_data_cleanup.m
% PURPOSE
%   Function to access a range of tools to clean-up timeseries data
% USAGE
%   ct_data_cleanup(mobj,src)
% INPUTS
%   muicat - handle to CoastalTools Dataset catalogue to allow access to data
%   src - menu object used to identify which function to be used
% OUTPUT
%   results are saved to the relevant CoastalTools class object
% NOTES
%   List of cleanup functions available:
%       Concatenate two timeseries
%       Resample timeseries
%       Delete multiple profiles
%       Edit or Delete profile in timeseries
%
% Author: Ian Townend
% CoastalSEA, (c) 2018, 2020 modified for use in MUI.CoastlTools
%--------------------------------------------------------------------------
%
    switch src.Text
        case 'Concatenate two timeseries'
            concatenate_ts(muicat);
        case 'Resample timeseries'
            resample_ts(muicat);
        case 'Patch timeseries'
            patch_ts(muicat);
        case 'Trim timeseries'
            trim_ts(muicat);
        case 'Delete multiple profiles'
            delete_profile_ts(muicat);
        case 'Edit or Delete profile in timeseries'
            edit_delete_profile(muicat);
    end
end
%%
function concatenate_ts(muicat)
    %concatenate two time series and write result as a new timeseries
    %get first time series
    seltype = [];
    promptxt = 'Select first timeseries';
    [caserec,ok] = selectTSdata(mobj,seltype,promptxt);
    if ok<1, return; end  %user cancelled    
    [tsc1,caseobj] = getCaseDataSet(mobj.Cases,mobj,caserec); 
    [classhandle,id_class,~] = getCaseDataID(mobj.Cases,[],caserec);   
    
    %get second time series
    promptxt = 'Select second timeseries';
    [caserec2,ok] = selectTSdata(mobj,seltype,promptxt);
    if ok<1, return; end  %user cancelled    
    [tsc2,~] = getCaseDataSet(mobj.Cases,mobj,caserec2);   
    
    %order in sequence and find any overlap    
    if datenum(tsc1.TimeInfo.StartDate)>datenum(tsc2.TimeInfo.StartDate)
        %tsc2 is before tsc1 so swap them around
        temp = tsc1;
        tsc1 = tsc2;
        tsc2 = temp;
        clear temp
    end
    %now get start and end time of both timeseries
    tsc1_startime = datenum(tsc1.TimeInfo.StartDate);
    tsc1_endtime = tsc1_startime+tsc1.TimeInfo.End;
    tsc2_startime = datenum(tsc2.TimeInfo.StartDate);        
    tsc2_endtime = tsc2_startime+tsc2.TimeInfo.End;
    
    if tsc1_endtime>tsc2_startime
        %there is an overlap
        
        
        quest = 'Use which timeseries for overlap?';
        tsn = questdlg(quest,'Concatenate timeseries',...
                                    tsc1.Name,tsc2.Name,tsc1.Name);
        promptxt = {'Adjust end of TS1', 'Adjust start of TS2'};
        defaults = {datestr(tsc1_endtime),datestr(tsc2_startime)};
        values = inputdlg(promptxt,'Option to adjust',1,defaults);
        
        if ~isempty(values)
            tsc1_endtime = datenum(values{1});
            tsc2_startime = datenum(values{2});
        end

        switch tsn
            case tsc1.Name
                %adjust tsc1 in case user has changed end date
                tsc1name = tsc1.Name;                
                tsc1 = getsampleusingtime(tsc1,datestr(tsc1_startime),...
                                                    datestr(tsc1_endtime));
                tsc1.Name = tsc1name;
                %trim tsc2 to date after tsc1_endtime
                tsc2name = tsc2.Name;
                offset = TSDataSet.ts_interval(tsc2,'First interval');
                %contrary to documentation, timeseries functions cannot handle datenums
                tsc1_endtime  = datestr(tsc1_endtime+offset);
                tsc2_endtime  = datestr(tsc2_endtime);
                tsc2 = getsampleusingtime(tsc2,tsc1_endtime,tsc2_endtime);
                tsc2.Name = tsc2name;
                switchtime = tsc1_endtime;                
            case tsc2.Name
                %adjust tsc2 in case user has changed start date
                tsc2name = tsc2.Name;
                tsc2 = getsampleusingtime(tsc2,datestr(tsc2_startime),...
                                                datestr(tsc2_endtime));
                tsc2.Name = tsc2name;
                %trim tsc1 to date before tsc2_startime
                tsc1name = tsc1.Name;
                offset = TSDataSet.ts_interval(tsc1,'First interval');                
                tsc1_startime = datestr(tsc1_startime);        
                tsc2_startime = datestr(tsc2_startime-offset);
                tsc1 = getsampleusingtime(tsc1,tsc1_startime,tsc2_startime);
                tsc1.Name = tsc1name;
                switchtime = tsc2_startime;                
        end
    else
        tsn = 'No overlap';
        switchtime = 'timeseries';
    end
    
    %concatenate single timeseries or two tscollections with same variables
    quest = 'Do you want to concatenate a single variable or all variables?';
    answer = questdlg(quest,'Concatenate timeseries',...
                                    'Single','All','Cancel','All');    
    switch answer
        case 'Cancel'
            return
        case 'Single'
            idx1 = 1; idx2 = 1;
            varnames = gettimeseriesnames(tsc1);
            if length(varnames)>1
                [idx1,ok] = listdlg('Name','TS1 options', ...
                    'PromptString','Select TS1 variable:', ...
                    'SelectionMode','single', ...
                    'ListString',varnames);
                if ok<1, idx1 = 1; end
            end
            ts1 = tsc1.(varnames(idx1));
            varnames = gettimeseriesnames(tsc2);
            if length(varnames)>1
                [idx2,ok] = listdlg('Name','TS2 options', ...
                    'PromptString','Select TS2 variable:', ...
                    'SelectionMode','single', ...
                    'ListString',varnames);
                if ok<1, idx2 = 1; end
            end
            ts2 = tsc2.(varnames(idx2));
            new_ts = append(ts1,ts2);   
            %add meta-data to timeseries
            new_ts.Name = varnames(idx2);
            new_ts.UserData = sprintf('Series used for overlap: %s',char(tsn));
            new_ts.DataInfo.UserData = sprintf('%s to %s',char(tsn),switchtime);
            new_tsc = tscollection(new_ts);
        case 'All'
            new_tsc = vertcat(tsc1,tsc2);
    end
    %save as a new record
    RecName = sprintf('%s and %s',tsc1.Name,tsc2.Name);            
    new_tsc.Name = RecName;
    new_tsc.TimeInfo.UserData = sprintf('%s to %s',char(tsn),switchtime);
    id_rec = length(caseobj.mtsc)+1;
    caseobj.mtsc{id_rec} = new_tsc;
    Results.saveResults(mobj,classhandle,caseobj,id_class);
     
end
%%
function resample_ts(muicat)
    %resample a timeseries or a tscollection to a different time interval. 
    %For timeseries nan values are kept. for a tscollection interolation 
    %is used to fill any gaps in each timeseries that make up the collection.
    %
    %select record to be used
    seltype = [];
    promptxt = 'Select a timeseries';
    [caserec,ok] = selectTSdata(mobj,seltype,promptxt);
    if ok<1, return; end  %user cancelled    
    [tsc1,caseobj] = getCaseDataSet(mobj.Cases,mobj,caserec); 
    [classhandle,id_class,~] = getCaseDataID(mobj.Cases,[],caserec); 
    
    %concatenate single timeseries or two tscollections with same variables
    quest = 'Do you want to resample a single variable or all variables?';
    answer = questdlg(quest,'Resample timeseries',...
                                    'Single','All','Cancel','All');   
    switch answer
        case 'Single'                           
            %select variable to be interpolated
            varnames = gettimeseriesnames(tsc1);
            if length(varnames)>1
                [idx,ok] = listdlg('Name','TS options', ...
                                    'PromptString','Select TS variable:', ...
                                    'SelectionMode','single', ...
                                    'ListString',varnames);
                if ok<1, return; end
            else
                idx = 1;
            end
            sts = tsc1.(varnames{idx});  %timeseries of selected variable
            %idnan = isnan(sts.Data);     %index for all nan values in source
            ptxt = varnames{idx};
        case 'All'
            sts = tsc1; 
            ptxt = 'All';
        otherwise
            return
    end
    tint = [];
    while isempty(tint)
        tint = get_timeinterval(sts);
    end
    %now resample the selected timeseries or tscollection
    startime = datetime(sts.TimeInfo.StartDate);
    endtime = datetime(startime+sts.TimeInfo.End);
    time = cellstr((startime:tint:endtime)');
    new_tsc = resample(sts,time); %applies to a timeseries or tscollection
    if isa(new_tsc,'timeseries')
        %new_tsc.Data(idnan) = NaN;       %restore nan values**************
        new_tsc = tscollection(new_tsc); %make a tscollection if timeseries
    end
    RecName = sprintf('%s resampled',tsc1.Name);           
    new_tsc.Name = RecName;
    new_tsc.TimeInfo.UserData = sprintf('%s, %s resampled at %s',tsc1.Name,ptxt,string(tint));

    id_rec = length(caseobj.mtsc)+1;
    caseobj.mtsc{id_rec} = new_tsc;
    Results.saveResults(mobj,classhandle,caseobj,id_class);
     
end
%%
function patch_ts(muicat)
    %replace NaN values in one timeseries with values from another ts
    %NB both timeseries must be the same length
    seltype = [];
    promptxt = 'Select primary timeseries';
    [caserec,ok] = selectTSdata(mobj,seltype,promptxt);
    if ok<1, return; end  %user cancelled    
    [tsc1,caseobj] = getCaseDataSet(mobj.Cases,mobj,caserec); 
    [classhandle,id_class,~] = getCaseDataID(mobj.Cases,[],caserec);   
    
    %get second time series
    promptxt = 'Select infill timeseries';
    [caserec2,ok] = selectTSdata(mobj,seltype,promptxt);
    if ok<1, return; end  %user cancelled    
    [tsc2,~] = getCaseDataSet(mobj.Cases,mobj,caserec2); 
    
    time = getabstime(tsc1);
    t1 = datenum(time);    
    t2 = datenum(getabstime(tsc2));

    idx1 = 1; idx2 = 1;
    varnames = gettimeseriesnames(tsc1);
    if length(varnames)>1
        [idx1,ok] = listdlg('Name','TS1 options', ...
            'PromptString','Select TS1 variable:', ...
            'SelectionMode','single', ...
            'ListString',varnames);
        if ok<1, idx1 = 1; end
    end
    ts1 = tsc1.(varnames(idx1));
    varnames = gettimeseriesnames(tsc2);
    if length(varnames)>1
        [idx2,ok] = listdlg('Name','TS2 options', ...
            'PromptString','Select TS2 variable:', ...
            'SelectionMode','single', ...
            'ListString',varnames);
        if ok<1, idx2 = 1; end
    end
    ts2 = tsc2.(varnames(idx2));
    
    data = ts1.Data;
    missing = isnan(data);
    [~,patch] = intersect(t2,t1(missing));
    data(missing)= ts2.Data(patch);
    
    new_ts = timeseries(data,time);  
    %add meta-data to timeseries
    new_ts.Name = ts1.Name;
    new_ts.UserData = sprintf('Series used for patch: %s',ts2.Name);
    new_tsc = tscollection(new_ts);
    
    %save as a new record
    RecName = sprintf('%s and %s',tsc1.Name,tsc2.Name);            
    new_tsc.Name = RecName;
    new_tsc.TimeInfo.UserData = sprintf('Series used for patch: %s',ts2.Name);
    id_rec = length(caseobj.mtsc)+1;
    caseobj.mtsc{id_rec} = new_tsc;
    Results.saveResults(mobj,classhandle,caseobj,id_class);
         
end
%%
function trim_ts(muicat)
    %allow user to adjust the start and end data of a timeseries
    %select record to be used
    seltype = [];
    promptxt = 'Select a timeseries';
    [caserec,ok] = selectTSdata(mobj,seltype,promptxt);
    if ok<1, return; end  %user cancelled    
    [tsc1,caseobj] = getCaseDataSet(mobj.Cases,mobj,caserec); 
    [classhandle,id_class,~] = getCaseDataID(mobj.Cases,[],caserec);
    
    %now get start and end time to use
    tsc1_startime = tsc1.TimeInfo.StartDate;
    tsc1_endtime = datestr(datenum(tsc1_startime)+tsc1.TimeInfo.End);
    promptxt = {'Adjust Start date', 'Adjust End date'};
    defaults = {datestr(tsc1_startime),datestr(tsc1_endtime)};
    values = inputdlg(promptxt,'Option to adjust',1,defaults);
    if ~isempty(values)
        tsc1_startime = values{1};
        tsc1_endtime = values{2};
    end
    
    tsc1name = tsc1.Name;                
    tsc1 = getsampleusingtime(tsc1,tsc1_startime,tsc1_endtime);
    tsc1.Name = tsc1name; 
    
    id_rec = length(caseobj.mtsc)+1;
    caseobj.mtsc{id_rec} = tsc1;
    Results.saveResults(mobj,classhandle,caseobj,id_class);
    
end
%%
function delete_profile_ts(muicat)
    %remove all profiles with fewer than N time steps available  
    classname = 'ctBeachProfileData';
    pobj = muicat.DataSets.(classname);
    
    prompt = {'Minimum number of time steps'};
    title = 'Delete short profile records';
    default = {num2str(0)};
    answer = inputdlg(prompt,title,1,default);
    nmin = str2double(answer{1});                
    %find all profiles of a selected format
    idfall = [pobj.idFormat]';
    ntype = length(unique(idfall));
    if ntype>1
        [idformat,ok] = listdlg('PromptString','Select a file format',...
                    'SelectionMode','single','ListSize',[300,100],...                    
                    'ListString',pobj(1).DataFormats(:,1));
        if ok<1, return; end %user cancelled
    else
        idformat = 1;
    end
    idfsel = logical(idfall==idformat);

    %find all records with less than tmin timesteps
    nrec = length(pobj);
    ntstep = zeros(nrec,1);
    for i=1:nrec
       ntstep(i,1) = height(pobj(i).Data.Dataset.DataTable);
    end

    idx = ntstep<nmin & idfsel;
    if any(idx)
        deleteID = [pobj(idx).CaseIndex];
        msgbox(sprintf('%d profiles deleted',length(deleteID)));
    else
        msgbox('No profiles below threshold set');
        return
    end

    scenariolist = muicat.Catalogue.CaseID;   
    %find the ids of the deletelist in the full scenariolist
    caserec = find(ismember(scenariolist,deleteID));
    deleteCases(muicat,caserec)  
end
%%
function edit_delete_profile(muicat)
    %edit or delete a single profile for a timeseries record
    isok = 1;    
    classname = 'ctBeachProfileData';
    datasetname = 'Dataset';
%     pobj = muicat.DataSets.(classname);
    
    promptxt = 'Select profile timeseries (Cancel to quit)';
    while isok>0
        [caserec,isok] = selectRecord(muicat,'PromptText',promptxt,...
                              'CaseClass',classname,'ListSize',[150,250]);
        if isok<1, return; end %user cancelled  

        [cobj,classrec,~] = getCase(muicat,caserec);
        dst = cobj.Data.(datasetname);

        ptime = dst.RowNames;
        hfig = figure('Name','Delete individual profiles','Tag','PlotFig',...
            'Units','normalized');
        hfig.Position(1) = 0.05; %move figure to top left
        hfig.Position(2) = 1-hfig.Position(4)-0.1;
        hax = axes(hfig); %#ok<LAXES>
        hold on
        for j=1:length(ptime)
            elevation = dst.Elevation(j,:);
            chainage = dst.Chainage(j,:);    
            pcolor = [0.9,0.9,0.9];
            plot(chainage,elevation,'Color',pcolor,'Tag','Background')                
        end
        ok = 1;
        select = questdlg('Select profiles or scroll all?',dst.Description,...
                                    'Select','Scroll','Cancel','Select');
        if strcmp(select,'Select')
            while ok>0
                [idx,ok] = listdlg('Name',dst.Description, ...
                    'PromptString','Select a time step:', ...
                    'SelectionMode','single', ...
                    'ListString',ptime);
                if ok>0
                    [dst,ptime] = plotdeleteprofile(dst,ptime,idx,hax);
                    muicat.DataSets.(classname)(classrec).Data.Dataset = dst;
                end
            end
        elseif strcmp(select,'Scroll')
            nrec = length(ptime); ndel = 0;
            for i=1:nrec
                %length of ptime changes if profiles deleted
                if length(ptime)<nrec
                   nrec = nrec-1;
                   ndel = ndel+1;
                end 
                ip = i-ndel;
                %
                if ip<=nrec
                    [dst,ptime,ok] = plotdeleteprofile(dst,ptime,ip,hax);
                    if ok>0
                        muicat.DataSets.(classname)(classrec).Data.Dataset = dst;
                    else
                        break
                    end
                end
            end
        end
        hold off
        close(hfig)
    end

        %nested function to plot a profile and allow user to delete it
        function [dst,ptime,ok] = plotdeleteprofile(dst,ptime,idx,hax)
            ok = 1;
            y = dst.Chainage(idx,:)'; %only one profile
            z = dst.Elevation(idx,:)'; 
            editline = findobj(hax,'Tag','EditLine');
            if ~isempty(editline)
                delete(editline)
            end
            plot(hax,y,z,'Tag','EditLine');
            legend(sprintf('Profile: %s, Date: %s',dst.Description,...
                                        datestr(ptime(idx),'dd-mmm-yyyy')));

            action = questdlg('Edit or Delete this profile?',dst.Description,...
                                    'Edit/Delete','Next','Cancel','Next');
            if strcmp(action,'Cancel')   
                ok = 0;
            elseif strcmp(action,'Edit/Delete')
                %option to edit or delete selected profile
                answer = questdlg('Edit or Delete this profile?',dst.Description,...
                                    'Edit','Delete','Cancel','Edit');
                if strcmp(answer,'Delete')
                    dst.DataTable(idx,:) = [];
                    ptime(idx) = [];
                elseif strcmp(answer,'Edit')
                    %zero chainage can be repeated in some profiles
                    [isall,idx] = isunique(y,false);
                    if ~isall                 %add a small offset to duplicate value
                        y(idx) = y(idx)+eps;  %only works if there is just one duplicate
                    end
                    idd = ~isnan(y);  %remove nans
                    %create tablefigure for editing data
                    oldtable = table(z(idd),'RowNames',string(y(idd)));
                    title = sprintf('Edit profile'); 
                    txt1 = 'To edit select a cell, amend value and press return, or select another cell';
                    header = sprintf('%s\nProfile: %s',txt1,dst.Description);    
                    but.Text = {'Save','Cancel'}; %labels for tab button definition
                    newtable = tablefigureUI(title,header,oldtable,true,but);
                    if isempty(newtable), return; end  %user cancelled  
                    
                    dst.Elevation(idx,:) = newtable{:,:};
                    nvar = length(dst.VariableQCflags);
                    dst.VariableQCflags(1:nvar) = repmat({'qc'},1,nvar);
                end
            end
        end
end

%%
function tint = get_timeinterval(ats)
    %prompt user for the time interval to be used and return as a duration
    %works if ats is a timeseries or a tscollection
    varunits = ats.TimeInfo.Units;
    tint1 = TSDataSet.ts_interval(ats,'First interval');
    stint1 = cellstr(tint1,varunits(1));
    stint2 = cellstr(TSDataSet.ts_interval(ats,'Mean'),varunits(1));
    stint3 = cellstr(TSDataSet.ts_interval(ats,'Mode'),varunits(1));
    p1 = sprintf('First interval %s, Mean %s, Mode %s',stint1{1},stint2{1},stint3{1});
    
    p2 = sprintf('Duration of time steps (%s)',varunits);
    prompt = {sprintf('%s\n%s',p1,p2),'Units (days, hours, minutes, seconds)'};
    title = 'Define interval';
    numlines = 1;
    default = {num2str(datenum(tint1),4),varunits};    
    answer = inputdlg(prompt,title,numlines,default);
    tint = str2double(answer{1}); %new interval
    switch answer{2}
        case 'days'
            tint = days(tint);
        case 'hours'
            tint = hours(tint);
        case 'minutes'
            tint = minutes(tint);    
        case 'seconds'
            tint = seconds(tint);
        otherwise
            msgbox('Units can only be "days","hours","minutes" or "seconds"');
            tint = [];
            return;
    end
end
% %%
% function [recnum,ok] = selectTSdata(mobj,seltype,promptxt)
%     %prompt user to select a timeseries data set
%     %prompt user to select an existing form model to modify
%     % seltype - handle of class to use in subselection
%     %        - or data type (model, data, etc)            
%     [recnum,~,~,ok] = ScenarioList(mobj.Cases,seltype,...
%                 'PromptText',promptxt,'ListSize',[300,100]);
% %     if ok<1
% %         warndlg('No selection made');
% %         return;
% %     end
% end
