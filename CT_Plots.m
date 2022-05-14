classdef CT_Plots < muiPlots
%
%-------class help---------------------------------------------------------
% NAME
%   CT_Plots.m
% PURPOSE
%   Class to implement the generation of a range of plot types in
%   CoastalTools
% NOTES
%   inherits from muiPlots
% SEE ALSO
%   called from CT_PlotsUI.m, which defines the selection and settings in
%   properties UIselection and UIsettings. Uses muiCatalogue to access
%   data. 
%
% Author: Ian Townend
% CoastalSEA (c)July 2021
%--------------------------------------------------------------------------
%
    properties (Transient)
        %inherits the following properties from muiPlots
        %Plot            %struct array for:
                         %   FigNum - index to figures created
                         %   CurrentFig - handle to current figure
                         %   Order - struct that defines variable order for
                         %           plot type options (selection held in Order)

        %ModelMovie      %store animation in case user wants to save
        %UIsel           %structure for the variable selection made in the UI
        %UIset           %structure for the plot settings made in the UI
        %Data            %data to use in plot (x,y,z)
        %TickLabels      %struct for XYZ tick labels
        %AxisLabels      %struct for XYZ axis labels
        %Legend          %Legend text
        %MetaData        %text summary of primary variable selection
        %Title           %Title text
        %Order           %order of variables for selected plot type
        %idxfig          %figure number of the current figure
    end
%%
    methods 
        function obj = CT_Plots
            %types of plot avaiable based on number of dimensions
            obj.Plot.FigNum = [];
            obj.Plot.CurrentFig = [];  
            obj.Plot.Order = CT_Plots.setVarOrder;
        end      
    end
%%    
    methods (Static)
        function getPlot(gobj,mobj)
            %get existing instance or create new class instance
            if isa(mobj.mUI.Plots,'CT_Plots') 
                obj = mobj.mUI.Plots;    %get existing instance          
                clearPreviousPlotData(obj);
            else
                obj = CT_Plots;                   %create new instance
            end

            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;
            
            %set the variable order for selected plot type
            obj.Order = obj.Plot.Order.(obj.UIset.callTab);
            
            %modify the legend text to include case and variable only
            legformat.idx = [1 0 1];
            legformat.text = {};
            
            dtype = 'array';
            switch obj.UIset.callTab
                case 'Timeseries'
                    addDefaultSelection(obj,'Time');
                case 'Profiles'     
                    addProfileChainage(obj,mobj);                    
                    dtype = 'dstable';
                    if obj.UIsel(2).property==1
                        obj.UIset.Type.String = 'contourf';
                        setMultiProfilePlot(obj,mobj)
                        return;
                    end
                case 'Rose'
                    addDefaultSelection(obj,'Dir',mobj);                   
                    legformat = repmat(legformat,1,2);
                    legformat(2).idx = [0 0 1]; %format as Case(Variable)Direction
                case {'2D','2DT'}
                    check4chainage(obj,mobj);   %check for use of chainage
                    legformat = repmat(legformat,1,2);
                case {'3D','3DT'}
                    check4chainage(obj,mobj);   %check for use of chainage
                    legformat = repmat(legformat,1,3);
            end

            %get the data to be used in the plot
            ok = getPlotData(obj,mobj.Cases,dtype,legformat);
            if ok<1, return; end %data not found
            isvalid = checkdimensions(obj);
            if ~isvalid, return; end
            
            if ~isempty(obj.UIset.Type) && strcmp(obj.UIset.Type.String,'User')
                user_plot(obj,mobj);  %pass control to user function
            else
                %generate the plot
                setPlot(obj,mobj);
            end
        end
     end
%%   
    methods (Access=protected)       
        function callPlotType(obj)
            %call the function specific to the selected plot type           
            switch obj.UIset.callTab        %call function based on Tab
                case {'Timeseries','Rose','2D'}                 
                    switch obj.UIset.callButton %and Tab Button used
                        case 'New'              %create new 2D plot
                            %check whether cartesian or polar plot
                            if obj.UIset.Polar
                                newPolarplot(obj);                              
                            else
                                new2Dplot(obj);  
                            end
                        case 'Add'              %add variable to 2D plot
                            if strcmp(obj.UIset.Type,'bar')
                                addBarplot(obj);
                            elseif obj.UIset.Polar
                                addPolarPlot(obj)  
                            else
                                add2Dplot(obj);
                            end
                        case 'Delete'           %delete variable from 2D plot
                            del2Dplot(obj);
                    end
                case 'Profiles'
                    setProfilePlot(obj)
                case '3D'
                    new3Dplot(obj);
                case '4D'
                    new4Dplot(obj);
                case {'2DT','3DT','4DT'}
                    newAnimation(obj);
                otherwise
                    warndlg('Could not find plot option in callPlotType');
            end
        end
%%
        function addDefaultSelection(obj,Prop,mobj)
            %add a selection that is included by default (eg time or
            %direction)
            obj.UIsel(2) = obj.UIsel(1);
            obj.UIsel(1).xyz = [true false];
            obj.UIsel(2).xyz = ~obj.UIsel(1).xyz; %invert logical array

            if strcmp(Prop,'Time')
                obj.UIsel(2).variable = 1;
                obj.UIsel(2).property = 2;
                trange = obj.UIsel(1).dims.value;
                obj.UIsel(2).range = trange;
                obj.UIsel(2).desc = 'Time';                
                obj.UIsel(2).dims = struct('name','','value',[]);
                obj.UIset.Polar = false;
            elseif strcmp(Prop,'Dir')
                crec = obj.UIsel(1).caserec;
                dset = obj.UIsel(1).dataset;
                dst = getDataset(mobj.Cases,crec,dset);
                idvar = find(contains(dst.VariableNames,'Dir'));
                obj.UIsel(2).variable = idvar;
                varange = dst.VariableRange.(dst.VariableNames{idvar});
                obj.UIsel(2).range = var2range(varange);                
                obj.UIsel(2).desc = 'Direction';
                obj.UIset.Polar = true;
            end
        end
%%
        function addProfileChainage(obj,mobj)
            %add the profile chainage variable to the selection 
            for i=1:length(obj.UIsel)
                if obj.UIsel(i).property==2
                    %selection is time so adjust selection accordingly
                    addDefaultSelection(obj,'Time',mobj);
                elseif obj.UIsel(i).caserec>0
                    %selection is a profile variable - add chainage
                    idvar = obj.UIsel(i).variable;
                    dst = getDataset(mobj.Cases,obj.UIsel(i).caserec,1);
                    fnames = dst.VariableNames;
                    idchain = find(strcmp(fnames,'Chainage'));
                    obj.UIsel(i).variable = [idvar,idchain];
                end
            end
        end
%%
        function check4chainage(obj,mobj)
            %if noDim is being used for Chainage, reassign selection 
            %NB - specific to ctBeachProfileData beach profiles
            crec = obj.UIsel(1).caserec;
            dset = obj.UIsel(1).dataset;
            if ~strcmp(mobj.Cases.Catalogue{crec,3},'ctBeachProfileData')
                return;
            end
            dst = getDataset(mobj.Cases,crec,dset);
            varnames = dst.VariableNames;
            for i=2:length(obj.UIsel)
                ui = obj.UIsel(i);
                if length(ui.dims)>1 && strcmp(ui.dims(2).name,'noDim1') && ...
                        ui.property==3
                    %if profile change second variable to chainage 
                    idx = find(strcmp(varnames,'Chainage')); %assumes Chainage is key word
                    obj.UIsel(i).variable = idx;
                    obj.UIsel(i).property = 1;
                    obj.UIsel(i).dims = obj.UIsel(1).dims;
                end
            end
        end
%%
%--------------------------------------------------------------------------
% Code to handle beach profiles
%--------------------------------------------------------------------------
        function setProfilePlot(obj)
            %generate new profile plot in figure            
            if isempty(obj.Data.X) 
                %plot multiple X-Z profiles
                getDistPlot(obj);  
            elseif isempty(obj.Data.X.DataTable)
                %user selected chainage dimension rather than time
                idx = obj.Plot.FigNum==obj.Plot.CurrentFig.Number;
                hf = findobj('Number',obj.Plot.FigNum(idx));
                delete(hf)
                obj.Plot.FigNum(idx)=[];
            elseif isdatetime(obj.Data.X.DataTable{1,1})
                %plot surface plot of time variation at one location
                timeProfilePlot(obj);
            else
                %two different profiles selected - plot multi-profiles
                %called directly in getPlot because needs mobj
            end
        end
%%
        function getDistPlot(obj)
            %definitions of layout for distance plot (for profle data)
            % obj - plot class instance
            % scr - button used to plot (New, Add, etc)
            %profile data assigned as x=chainage, y=elevation, z=time
            ok = 1;
            idx = obj.Plot.FigNum==obj.idxfig;
            hfig = findobj('Number',obj.Plot.FigNum(idx));
            xall = obj.Data.Y.DataTable{:,2};    %chainage in extracted table
            yall = obj.Data.Y.DataTable{:,1};    %variable in extracted table
            tall = obj.Data.Y.RowNames; %time   
            talltxt = cellstr(tall);
            %get the legend based on user selection (New, Add)
            tlist = getSelectionList(obj,talltxt);
            legtxt = getlegendText(obj,hfig,talltxt);
            %
            while ok>0
                [h_dlg,ok] = listdlg('Name','Plot profile', ...
                    'PromptString','Select time','ListSize',[180,300],... ...
                    'SelectionMode','single','ListString',tlist);
                if ok==0, return; end
                selection = tlist(h_dlg,:);
                if strcmp(selection,'<   All profiles   >') || ...
                    strcmp(selection,'<  Winter profiles >') || ...
                    strcmp(selection,'<  Summer profiles >') 
                
                    %plot all profiles for a single location
                    ok = plotAllProfiles(obj,hfig,xall,yall,tall,selection);
                else
                    %add and delete selected profiles individually 
                    tlist{h_dlg} = sprintf('%s - plotted',tlist{h_dlg});
                    plotSingleProfile(obj,hfig,...
                                             xall,yall,legtxt,h_dlg);
                end
            end
        end   
%%
        function tlist = getSelectionList(obj,tall)
            %get the list of options depending on the call (New, add)
            tlist = tall;            %initialise date selection list
            %
            if strcmp(obj.UIset.callButton,'New')
                adtxt = {'<   All profiles   >';'<  Winter profiles >';...
                                                '<  Summer profiles >'};
                tlist = vertcat(tall,adtxt);
            end
        end
%%
        function legtxt = getlegendText(obj,hfig,tall)
            %set the legend and update existing plot if different profile added
            datetxt = split(tall); 
            profname = obj.Data.Y.Description;
            if strcmp(obj.UIset.callButton,'New')
                obj.AxisLabels.X = 'Chainage (m)';     
                obj.Title = sprintf('Profile: %s',profname);
                hfig.UserData = profname; 
                legtxt = datetxt(:,1);
            elseif strcmp(obj.UIset.callButton,'Add') && ~strcmp(profname,hfig.UserData)
                %first addition of a different profile
                hplots = findobj(hfig,'Type',obj.UIsel.PlotType);
                if isempty(hplots), return; end
                pname = hfig.UserData;
                for j=1:length(hplots)
                    pdate = hplots(j).DisplayName;
                    if ~strcmp(pdate,'0 mOD')
                        hplots(j).DisplayName = sprintf('%s: %s',pname,pdate);
                    end
                end
                figax = hfig.CurrentAxes; 
                figax.Title.String = 'Multiple profile plot';
                for i=1:length(tall)
                    legtxt{i} = sprintf('%s: %s',profname,datetxt{i,1}); 
                end
                hfig.UserData = 'AddMore';
            elseif strcmp(obj.UIset.callButton,'Add') &&  strcmp(hfig.UserData,'AddMore')
                %further additions of more than one profile location
                for i=1:length(tall)
                    legtxt{i} = sprintf('%s: %s',profname,datetxt{i,1}); 
                end
            else
                %add to the same profile used for New
                legtxt = datetxt(:,1);
            end
        end
%%
        function [subx,suby,subt] = getSubsetProfiles(~,x,y,dt,season)
            %subsample beach profiles and return winter or summer profiles
            %winter is from Oct-Mar and Summer is Apr-Sep
            switch season
                case 'winter'
                   idt = month(dt)<4 | month(dt)>9;            
                case 'summer'
                   idt = month(dt)>3 & month(dt)<10;           
            end
            subt = dt(idt); subx = x(idt,:); suby = y(idt,:);
        end
%%
        function ok = plotAllProfiles(obj,hfig,xall,yall,tall,selection)
            %plot all profiles for a single location
            nrows = size(xall,1);
            %give user option to omit legend if lots of profiles
            excleg = false;
            if nrows>10 
                answer = questdlg('Include legend (>10):','Long legend',...
                            'Yes','No','Yes');
                if strcmp(answer,'No')                            
                    excleg = true;  
                end
            end
            %subsample the data set if winter or summer selected
            if strcmp(selection,'<  Winter profiles >')
                [x,y,tall] = getSubsetProfiles(obj,xall,yall,...
                                                    tall,'winter');
                 obj.Title = [obj.Title,' (Winter)'];                               
            elseif strcmp(selection,'<  Summer profiles >')
                [x,y,tall] = getSubsetProfiles(obj,xall,yall,...
                                                tall,'summer');
                 obj.Title = [obj.Title,' (Summer)'];                           
            else
                x = xall; y = yall;
            end
            %plot the data as a seris of lines using newXYplot and addXYplot
            nline = size(x,1);
            set(groot,'defaultAxesColorOrder',jet(nline))
            ok = 0;
            obj.Data.X = x(1,:);
            obj.Data.Y = y(1,:);
            obj.Data.Z = [];
            datetxt = split(cellstr(tall));
            obj.Legend = datetxt{1,1};
            new2Dplot(obj)
            for ip=2:nline
                obj.Data.X = x(ip,:);
                obj.Data.Y = y(ip,:);
                obj.Legend = datetxt{ip,1};
                add2Dplot(obj)                          
            end
            set(groot,'defaultAxesColorOrder','remove')
            if excleg
                %remove the legend if user requested
                hl = findobj(hfig,'Type','legend');
                delete(hl); 
                colormap jet;
                colorbar('Ticks',[0,1],...
                        'TickLabels',{datetxt{1,1},datetxt{end,1}});
            end
            figax = hfig.CurrentAxes; 
            hold(figax,'on')    
            addProfilePlotFeatures(obj,figax)
            hold(figax,'off')
        end
%%
        function plotSingleProfile(obj,hfig,xall,yall,legtxt,h_dlg)                                             
            %select profiles individually                  
            obj.Data.X = xall(h_dlg,:);
            obj.Data.Y = yall(h_dlg,:);
            obj.Data.Z = [];           
            obj.Legend = legtxt{h_dlg};
            switch obj.UIset.callButton
                case 'New'
                    new2Dplot(obj);
                    obj.UIset.callButton = 'Add';
                    figax = hfig.CurrentAxes; 
                    hold(figax,'on')    
                    addProfilePlotFeatures(obj,figax)
                    hold(figax,'off')
                case 'Add'
                    add2Dplot(obj);                         
                case 'Delete'
                    del2Dplot(obj);
            end
        end
%%
        function addProfilePlotFeatures(obj,figax)
            %add zero line to elevation and change y-axis labels for features
            if strcmp(obj.AxisLabels.Y,'FeatureCode')
                %change Y-labels for feature code data
                SFC = SFcodes();
                catnames = SFC(:,2);
                figax.YTick = 1:length(catnames);
                figax.YTickLabel = catnames;
            elseif strcmp(obj.AxisLabels.Y,'Elevation (mOD)')
                %add OmOD for elevation data
                plot(figax,xlim,0*[1 1],'-.k','LineWidth',0.2,...
                            'DisplayName','0 mOD');
            end
        end
%%
        function timeProfilePlot(obj)
            %plot surface plot of time variation at one location
            idchain = strcmp(obj.Data.Y.DataTable.Properties.VariableNames,...
                                'Chainage');
            xall = obj.Data.Y.DataTable{:,idchain}; %chainage in extracted table            
            xall(isnan(xall)) = 0;  %X and Y co-ordinates cannot be NaN
            dates = datenum(obj.Data.Y.RowNames);       %time
%             dates.Format = 'dd-MMM-yyyy';
            yall = repmat(dates,1,size(xall,2));
            zall = obj.Data.Y.DataTable{:,~idchain}; %variable in extracted table           
            xint = size(xall,2)-1;
            yint = length(dates)-1;
            Xtext = 'Chainage (m)';
            Ytext = 'Time';
            Title = sprintf('Profile: %s',obj.Data.Y.Description);
            ptype = 'contourf';
            muiPlots.xySurface(xall,yall,zall,xint,yint,Xtext,Ytext,...
                                     obj.Legend,Title,ptype);
            ax = gca;
            ax.YTickLabel = cellstr(datetime(ax.YTick,'ConvertFrom','datenum'));
        end
%%
        function setMultiProfilePlot(obj,mobj)
            %plot multi-profile. call bypassess muiPlots.setPlot so need to
            %get figure first and 
            %get the data to be used in the plot
            ok = getPlotData(obj,mobj.Cases,'dstable');
            if ok<1, return; end %data not found
            isvalid = checkdimensions(obj);
            if ~isvalid, return; end
            getFigure(obj); 
            %call the specific plot type requested
            multiProfilePlot(obj,mobj)
            %assign muiPlots instance to handle
            mobj.mUI.Plots = obj;
        end  
%%            
        function multiProfilePlot(obj,mobj)
            %plot multi-profile. 
            idx = obj.Plot.FigNum==obj.idxfig;
            hfig = findobj('Number',obj.Plot.FigNum(idx));            
            muicat = mobj.Cases.Catalogue;
            prec = find(strcmp(muicat.CaseClass,'ctBeachProfileData'));
            profnames = muicat.CaseDescription(prec);
            %check order of profiles along shoreline
            [~,~,idd] = shore_profile_order(mobj,prec);
            sortedprofnames = profnames(idd);
            startprof = obj.Data.Y.Description;
            endprof = obj.Data.X.Description;
            idst = find(strcmp(sortedprofnames,startprof));
            idnd = find(strcmp(sortedprofnames,endprof));
            %casrecs of ordered selection
            if idst>idnd
                temp = idst;
                idst = idnd; idnd = temp;
            end
            profnames = sortedprofnames(idst:idnd);
            profrec = prec(idd(idst:idnd));
            
            nprof= length(profrec);
            nrows = zeros(nprof,1); mcols = nrows;
            obj.Data = [];
            for i=1:nprof
                %Construct valid MATLAB identifiers from input strings
                fname = matlab.lang.makeValidName(profnames{i},'ReplacementStyle','delete');
                profnames{i} = fname;
                obj.Data.(fname) = getDataset(mobj.Cases,profrec(i),1);
                [nrows(i),mcols(i)] = size(obj.Data.(fname).Chainage);
            end
            mxch = max(mcols);                   %maximum length of chainage record
            variable = obj.UIsel(1).variable(1); %assume first selection defines variable to use
            times = get_profile_times(mobj,profrec);
            times.Format = 'dd-MMM-yyyy';
            
            %extractt he data from each profile dataset
            mxtm = length(times);             
            x = zeros(nprof,mxch,mxtm); 
            y = repmat(1:nprof,length(times),1,mxch);
            z = NaN(nprof,mxch,mxtm);
            for i=1:nprof
                    dst = obj.Data.(profnames{i});
                    var = dst.DataTable{:,variable};
                    ptime = obj.Data.(profnames{i}).RowNames;%time in first profile
                    %size of x and z not same for all profiles
                    zr = 1:mcols(i);
                    tr = find(ismember(times,ptime)); %ids of ptime in times    
                    z(i,zr,tr) = var'; %(y,t,z)
                    Ch = dst.Chainage;
                    Ch(isnan(Ch)) = 0;  %X and Y co-ordinates cannot be NaN
                    x(i,zr,tr) = Ch';
            end
            
            %reverse order of profiles if required
            ylabel = profnames;
            qstr = sprintf('Profiles are from %s to %s',...
                profnames{1},profnames{end});
            revbut = questdlg(qstr,'Reverese order','Use','Reverse','Use');
            if strcmp(revbut,'Reverse')
                x = flipud(x);
                ylabel = flipud(ylabel);
                z = flipud(z);
            end 

            tlist = cellstr(times,'dd-MMM-yyyy');
            %now re-order so that time is the first dimension - for animation
            x = shiftdim(x,2);
            z = shiftdim(z,2);
            
            %set up Title and Legend
            vardesc = dst.VariableDescriptions{variable};
            obj.Title = sprintf('%s from %s to %s',vardesc,startprof,endprof);
            obj.Legend = dst.VariableLabels{variable};
            
            tbut = questdlg(qstr,'Type of output','Plot @t','Animation','Plot @t');
            if strcmp(tbut,'Plot @t')
                isok = 1;
                while isok>0
                    %select individual survey times for plotting
                    if isok>1
                        %get an existing figure of create a new one
                        getFigure(obj);
                        idx = obj.Plot.FigNum==obj.idxfig;
                        hfig = findobj('Number',obj.Plot.FigNum(idx));
                    end
                    plotSnapShot(x,y,z);
                    isok = isok+1;
                end
            else                
                %create an animated sequence of surface plots
                plotAnimation(obj,x,y,z)
            end  
            
            %% nested function---------------------------------------------
            function plotSnapShot(x,y,z)
                %select individual survey times for plotting
                [idt,ok] = listdlg('Name','Plot profile', ...
                        'PromptString','Select time','ListSize',[180,300],... ...
                        'SelectionMode','single','ListString',tlist);
                if ok==0, isok = -1; delete(hfig); return; end

                xt = squeeze(x(idt,:,:));
                yt = squeeze(y(idt,:,:));
                zt = squeeze(z(idt,:,:));

                Xtext = 'Chainage (m)';
                Ytext = 'Alongshore Profile Location';
                endidx = length(obj.Title);
                if contains(obj.Title,' on ')
                    endidx = regexp(obj.Title,' on ')-1;
                end
                obj.Title = sprintf('%s on %s',obj.Title(1:endidx),tlist{idt});
                muiPlots.xySurface(xt,yt,zt,mxtm-1,mxch-1,Xtext,Ytext,...
                            obj.Legend,obj.Title,obj.UIset.Type.String);                             
                figax = hfig.CurrentAxes; 
                hold(figax,'on')   
                figax.YTick = 1:nprof;
                figax.YTickLabel = ylabel;
                hold(figax,'off')
                if isvalid(obj.Plot.CurrentFig)
                    obj.Plot.CurrentFig.Visible = 'on';
                end
            end
            %% nested function---------------------------------------------
            function plotAnimation(obj,x,y,z)
                %create an animated sequence of surface plots
                xint = nprof-1; yint = mxch-1;
                for idt=1:length(tlist)                    
                    xt = squeeze(x(idt,:,:));
                    yt = squeeze(y(idt,:,:));
                    zt = squeeze(z(idt,:,:));
                    [xq,yq,zq(idt,:,:)] = gridAnimationData(obj,xt,yt,zt,xint,yint);
                end
                
                obj.Data = [];
                obj.Data.X = xq(1,:);
                obj.Data.Y = yq(:,1);
                obj.Data.Z = zq;
                obj.Data.T = times;
                obj.UIset.Type.String = 'contourf';
                obj.AxisLabels.X = 'Chainage (m)';
                obj.AxisLabels.Y = 'Alongshore Profile Location';

                %modify the tick labels on the y-axis   
                obj.TickLabels.YTick =1:nprof;
                obj.TickLabels.YTickLabel = ylabel;
                obj.UIset.callTab = '3DT'; %use default tab name in muiPlots
                
                %create animation
                newAnimation(obj);
            end
        end
%%
        function [xq,yq,zq] = gridAnimationData(~,x,y,z,xint,yint)
            %regrid the data based on user selected interval
            wid = 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId';
            minX = min(min(x)); maxX = max(max(x));
            minY = min(min(y)); maxY = max(max(y));
            xint = (minX:(maxX-minX)/xint:maxX);
            yint = (minY:(maxY-minY)/yint:maxY);
            [xq,yq] = meshgrid(xint,yint);
            warning('off',wid)
             zq = griddata(x,y,z,xq,yq);
            warning('on',wid)
        end
    end
%%
%--------------------------------------------------------------------------
% Static CT_Plots functions
%--------------------------------------------------------------------------
    methods(Static, Access=protected)
        function varorder = setVarOrder()
            %struct that holds the order of the variables for different
            %plot types        
            varnames = {'Timeseries','Profiles','Rose','2D','3D','2DT','3DT',};
            %types of plot in 2,3 and 4D              
            d2 = {'Y','X'}; ts = d2; pr = ts; ro = ts;
            d3 = {'Z','X','Y'};
            %types of animaton in 2,3 and 4D        
            t2 = {'Y','T','X'};
            t3 = {'Z','T','X','Y'};
            varorder = table(ts,pr,ro,d2,d3,t2,t3,'VariableNames',varnames);
        end
    end
end