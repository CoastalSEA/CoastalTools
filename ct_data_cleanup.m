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
%       Concatenate two timeseries datasets
%       Resample timeseries dataset
%       Patch NaNs or gaps in one timeseries with data from another
%       Trim the ends of a timeseries data set
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
    promptxt = 'Select first timeseries';
    [caserec1,isok] = selectRecord(muicat,'PromptText',promptxt,...
                                                    'ListSize',[150,250]);
    if isok<1, return; end %user cancelled 
    dst1 = copy(getDataset(muicat,caserec1,1));
    classname = muicat.Catalogue.CaseClass(caserec1);
    range1 = dst1.RowRange;
    [varname1,vidx1] = getVariable(dst1);
    
    %get second time series
    promptxt = 'Select second timeseries';    
    [caserec2,isok] = selectRecord(muicat,'PromptText',promptxt,...
                                                    'ListSize',[150,250]);
    if isok<1, return; end %user cancelled  
    dst2 = copy(getDataset(muicat,caserec2,1));
    range2 = dst2.RowRange;
    [varname2,vidx2] = getVariable(dst2);
    
    %order in sequence and find any overlap    
    if range1{1}>range2{1}
        %dst2 is before dst1 so swap them around
        temp1 = range1;     temp2 = dst1;
        range1 = range2;    dst1 = dst2;
        range2 = temp1;     dst2 = temp2;
        clear temp1 temp2
    end
    
    switchtime = 'end of first ts';
    if range1{2}>range2{1}
        %there is an overlap
        promptxt = {'Adjust end of TS1', 'Adjust start of TS2'};
        defaults = {datestr(range1{2}),datestr(range2{1})};
        values = inputdlg(promptxt,'Option to adjust',1,defaults);
        if ~isempty(values)
            range1{2} = datetime(values{1});
            range2{1} = datetime(values{2});
        end
        %adjust dst1 to account far any change in end date
        startime = range1{1}-minutes(1);  %offset ensures selected 
        endtime = range1{2}+minutes(1);   %range is extracted
        if ~checkdates(startime,endtime,dst.RowRange), return; end
        
        timeidx = isbetween(dst1.RowNames,startime,endtime);
        dst1 = getDSTable(dst1,timeidx);        
        %trim dst2 to date after dst1 endtime
        timeidx = isbetween(dst2.RowNames,endtime,range2{2}+minutes(1));
        dst2 = getDSTable(dst2,timeidx);        
        switchtime = datestr(range1{2}); 
    end
    
    if width(dst1)>1
        dst1 = removevars(dst1,dst1.VariableNames(~vidx1));
    end
    
    if width(dst2)>1
        dst2 = removevars(dst2,dst2.VariableNames(~vidx2));
    end
    
    if ~strcmp(varname1,varname2) %check that variable names are the same
        quest = 'Which variable name do you want to use?';
        answer = questdlg(quest,'Select variable name',...
                                    varname1,varname2,varname1);
        if strcmp(answer,varname2)
            %vercat uses first dstable for dsproperties so amend to dst2
            dst1.VariableNames{1} = varname2;
            dst1.VariableDescriptions{1} = dst2.VariableDescriptions{1};
            dst1.VariableLabels{1} = dst2.VariableLabels{1};
        else
            dst2.VariableNames{1} = varname1;
        end
    end
    
    %concatenate a new dstable
    newdst = vertcat(dst1,dst2);
    if isempty(newdst), return; end %failed to concatenate dstables
    
    newdst.Description = sprintf('%s and %s',dst1.Description,dst2.Description);
    newdst.Source = sprintf('%s, switched at %s',newdst.Description,switchtime);
    %save results as a new Record in Catalogue
    type = convertStringsToChars(muicat.Catalogue.CaseType(caserec1));    
    heq = str2func(classname);
    obj = heq();  %new instance of class object
    obj.Data.Dataset = newdst;   
    setCase(muicat,obj,type);
    getdialog(sprintf('Concatenated %s',newdst.Description));
end
%%
function resample_ts(muicat)
    %resample a timeseries to a different time interval. 
    datasetname = 'Dataset';   %uses default dataset name
    %select record to be used    
    promptxt = 'Select profile timeseries (Cancel to quit)';
    [caserec,isok] = selectRecord(muicat,'PromptText',promptxt,...
                                                    'ListSize',[150,250]);
    if isok<1, return; end %user cancelled  
    [cobj,~,catrec] = getCase(muicat,caserec);
    dst = cobj.Data.(datasetname);
    
    %get the old and new times for resampling
    tint = [];
    while isempty(tint)
        tint = get_timeinterval(dst,0);
    end
    oldtime = dst.RowNames;
    stend = dst.RowRange;
    newtime = (stend{1}:tint:stend{2})';
    
    %get the variables to be resampled
    varnames = dst.VariableNames;
    quest = 'Do you want to resample a single variable or all variables?';
    answer = questdlg(quest,'Resample timeseries',...
                                    'Single','All','Cancel','Single');   
    switch answer
        case 'Single'                           
            %select variable to be interpolated            
            if length(varnames)>1
                [idx,ok] = listdlg('Name','TS options', ...
                                    'PromptString','Select TS variable:', ...
                                    'SelectionMode','single', ...
                                    'ListString',varnames);
                if ok<1, return; end
            else
                idx = 1;
            end
            %now resample the selected timeseries
            ptxt = varnames{idx};
            dst = getDSTable(dst,[],idx);
            newvar = {interp1(oldtime,dst.(varnames{idx}),newtime,'linear','extrap')};            
        case 'All'
            ptxt = 'All';
            %now resample the selected timeseries
            nrec = length(varnames);
            newvar{1,nrec} = [];
            for i=1:nrec                
                if isinteger(dst.(varnames{i}))
                    oldvar = single(dst.(varnames{i}));
                else
                    oldvar = dst.(varnames{i});
                end
                newvar{i} = interp1(oldtime,oldvar,newtime,'linear','extrap');
            end
        otherwise
            return
    end
    
    %now assign new dataset to a dstable
    dsp = dst.DSproperties;
    newdst = dstable(newvar{:},'RowNames',newtime,'DSproperties',dsp);
    newdst.Description = sprintf('%s resampled',dst.Description);           
    newdst.Source = sprintf('%s, %s resampled at %s',dst.Description,ptxt,string(tint));
    
    %save results as a new Record in Catalogue
    type = convertStringsToChars(catrec.CaseType);
    heq = str2func(catrec.CaseClass);
    obj = heq();  %new instance of class object
    obj.Data.Dataset = newdst;   
    setCase(muicat,obj,type);
    getdialog(sprintf('Data resampled for: %s',catrec.CaseDescription));
end
%%
function patch_ts(muicat)
    %replace NaN values in one timeseries with values from another ts
    
    %get the first time series                                
    promptxt1 = 'Select primary timeseries';   
    [caserec1,isok] = selectRecord(muicat,'PromptText',promptxt1,...
                                                        'ListSize',[150,250]);
    if isok<1, return; end %user cancelled  
    dst1 = getDataset(muicat,caserec1,1);
    [tint,dst1] = get_timeinterval(dst1,1);  %if resampled then dst1 is a new dstable
    if isempty(tint), return; end
    t1 = dst1.RowNames;
    [varname1,vidx] = getVariable(dst1);
    ds1 = dst1.(varname1);

    stend = dst1.RowRange;
    newtime = (stend{1}:tint:stend{2})';
    %check whether the primary data set has time defined with NaNs or has
    %missing time rows
    if length(newtime)>length(t1)
        [~,patch] = intersect(newtime,t1);
        newvar = NaN(length(newtime),1);
        newvar(patch) = ds1;
    elseif length(newtime)<length(t1)
        warndlg('Case not handled in ct_data_cleanup.patch_ts')
        return;
    else
        newtime = t1; newvar = ds1;
    end
    
    %get the second time series
    promptxt2 = 'Select infill timeseries';
    [caserec2,isok] = selectRecord(muicat,'PromptText',promptxt2,...
                                                        'ListSize',[150,250]);
    if isok<1, return; end %user cancelled  
    dst2 = getDataset(muicat,caserec2,1);
    t2 = dst2.RowNames;
    varname2 = getVariable(dst2);
    ds2 = dst2.(varname2);
    
    %add the patch
    missing = isnan(newvar);
    patch = interp1(t2,ds2,newtime,'linear');
    newvar(missing)= patch(missing);
    
    %assign new dataset to a dstable
    dsp = dst1.DSproperties;
    dsp.Variables(~vidx) = [];
    newdst = dstable(newvar,'RowNames',newtime,'DSproperties',dsp);
    newdst.Description = sprintf('%s and %s',dst1.Description,dst2.Description);
    if iscell(dst1.Source)
        srctxt = dst1.Source{1};
    else
        srctxt = dst1.Source;
    end
    newdst.Source = {sprintf('%s patched with %s',srctxt,newdst.Description)};
    %save results as a new Record in Catalogue
    classname = muicat.Catalogue.CaseClass(caserec1);
    type = convertStringsToChars(muicat.Catalogue.CaseType(caserec1));
    heq = str2func(classname);
    obj = heq();  %new instance of class object
    obj.Data.Dataset = newdst;  
    setCase(muicat,obj,type);
    getdialog(sprintf('Patched %s',newdst.Description));         
end
%%
function trim_ts(muicat)
    %allow user to adjust the start and end data of a timeseries
    %select record to be used
    datasetname = 'Dataset';
    promptxt = 'Select profile timeseries (Cancel to quit)';
    [caserec,isok] = selectRecord(muicat,'PromptText',promptxt,...
                                                    'ListSize',[150,250]);
    if isok<1, return; end %user cancelled  
    [cobj,classrec,catrec] = getCase(muicat,caserec); %use getCase because need classrec
    classname = catrec.CaseClass; 
    dst = cobj.Data.(datasetname);
    
    %get start and end time to use and extract dataset
    values = editrange_ui(dst.RowRange);
    startime = datetime(values{1})-minutes(1);  %offset ensures selected 
    endtime = datetime(values{2})+minutes(1);   %range is extracted
    if ~checkdates(startime,endtime,dst.RowRange), return; end
  
    timeidx = isbetween(dst.RowNames,startime,endtime);    
    newdst = getDSTable(dst,timeidx);
    
    answer = questdlg('Update or Add new record?','Trim TS','Update','Add','Update');
    if strcmp(answer,'Update')
        muicat.DataSets.(classname)(classrec).Data.Dataset = newdst;
    else
        %save results as a new Record in Catalogue
        type = convertStringsToChars(catrec.CaseType);
        heq = str2func(catrec.CaseClass);
        obj = heq();  %new instance of class object
        obj.Data.Dataset = newdst;    
        setCase(muicat,obj,type);
        getdialog(sprintf('Data resampled for: %s',catrec.CaseDescription));
    end
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
    promptxt = 'Select profile timeseries (Cancel to quit)';
    while isok>0
        [caserec,isok] = selectRecord(muicat,'PromptText',promptxt,...
                              'CaseClass',{classname},'ListSize',[150,250]);
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
        %------------------------------------------------------------------
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
                    [isall,idy] = isunique(y,false);
                    if ~isall                 %add a small offset to duplicate value
                        y(idy) = y(idy)+eps;  %only works if there is just one duplicate
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
                    
                    dst.Elevation(idx,idd) = newtable{:,:};
                    nvar = length(dst.VariableQCflags);
                    dst.VariableQCflags(1:nvar) = repmat({'qc'},1,nvar);
                end
            end
        end
end
%%
function [tint,newdst] = get_timeinterval(dst,isresample)
    %prompt user for the time interval to be used and return as a duration
    %isresample -flag true if dst can be resampled; false when called by
    %resample. newvar returns input dst unless dst is resampled
    tint = []; newdst = dst;
    time = dst.RowNames;
    rowunit = dst.RowUnit;
    tint1 = ts_interval(time,rowunit,'First interval'); %returns a duration
    if contains(tint1.Format,':')
        rowunit = tint1.Format;
    end
    stint1 = cellstr(tint1,rowunit);
    
    tstep = unique(diff(time));   %unique timesteps as time durations
    if length(tstep)>1 && isresample
        %multiple timesteps in timeseries. resample to a single timestep
        action = questdlg('Multiple timesteps?','Time interval',...
                                    'Resample','Cancel','Resample');
        if strcmp(action,'Cancel'), return; end
        p1 = compose('%s; ',tstep)';
        p1 = sprintf('Time steps in selected timeseries\n%s',[p1{:}]);
        default = [char(tstep(1)),{rowunit}]; 
    else       
        stint2 = cellstr(ts_interval(time,rowunit,'Mean'),rowunit);
        stint3 = cellstr(ts_interval(time,rowunit,'Mode'),rowunit);
        p1 = sprintf('First interval %s, Mean %s, Mode %s',stint1{1},stint2{1},stint3{1});     
        default = {stint1{1},rowunit}; 
    end
    p2 = sprintf('Duration of time steps (%s)',rowunit);
    prompt = {sprintf('%s\n%s',p1,p2),'Units (days, hours, minutes, seconds)'};  
    title = 'Confirm interval';
    numlines = 1;       
    answer = inputdlg(prompt,title,numlines,default);
    if isempty(answer), return; end           %user cancelled
    tint = str2duration(answer{1},answer{2}); %new interval
    
    %if timeseries is to be resampled
    if (length(tstep)>1 && isresample) || (tint1~=tint)
        stend = dst.RowRange;
        newtime = (stend{1}:tint:stend{2})'; 
        varnames = dst.VariableNames;
        for i=1:length(varnames)
            if isinteger(dst.(varnames{i}))
                oldvar = single(dst.(varnames{i}));
            else
                oldvar = dst.(varnames{i});
            end
            newvar{i} = interp1(time,oldvar,newtime,'linear','extrap');
        end
        newdst = dstable(newvar{:},'RowNames',newtime,'DSproperties',dst.DSproperties);
        if iscell(dst.Source)    %restore Source desctiption
            newdst.Source = dst.Source;
        else
            newdst.Source = sprintf('Resampled: %s',dst.Source); 
        end
    end    
end
%%
function [varname,vidx] = getVariable(dst)
    %prompt user to select a variables from the dataset
    vars = dst.VariableNames;
    nrec = length(vars);
    idx = 1;
    if nrec>1
        [idx,ok] = listdlg('Name','TS options', ...
            'PromptString','Select TS variable:', ...
            'SelectionMode','single', ...
            'ListString',vars);
        if ok<1, idx = 1; end
    end
    varname = vars{idx};
    vidx = (1:nrec)==idx;
end
%%
function isok = checkdates(startime,endtime,RowRange)
    %check that start and end times not empty, are within range and in
    %correct order
    if isempty(startime) || isempty(endtime)
        msg = 'Incorrect date selection. No date selected for start or end';
        warndlg(msg);
        isok = false;
    else
        isok  = isvalidrange({startime,endtime},RowRange);
    end
end