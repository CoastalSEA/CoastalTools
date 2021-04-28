classdef Sim_BMV < muiDataSet
%
%-------class help---------------------------------------------------------
% NAME
%   Sim_BMV.m
% PURPOSE
%   Class description - Class for Sim_BMV model to run as a CoastalTools app
% USAGE
%   obj = Sim_BMV.runSimulation(mobj); %mobj is a handle to Main UI
%   obj = Sim_BMV.runForecast(mobj);
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:  
        TabOutput    %output struct with results from model. Struct fields:
                     % table - values saved for display on tabs
                     % headtxt - source and metatext to be used as table header
                     %BVI models also have a bvitable field
    end
    
    properties (Transient)
        UIsel           %structure for the variable selection made in the UI
        UIset           %structure for the plot settings made in the UI
    end
    
    methods (Access = private)
        function obj = Sim_BMV()                    
            %class constructor
        end
    end      
%%
    methods (Static)                
        function obj = runDataFitting(gobj,mobj)
            %function to run a simple 2D diffusion model
            obj = Sim_BMV;                           
            %now check that the input data has been entered
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); 
            %--------------------------------------------------------------------------
            % Model code 
            %--------------------------------------------------------------------------
            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;

            idx = [14,15,16,18];         %subset of variables to be checked
            [bvobj,prdst,wvdst,bmv,prwvs,bvrec] = getSimData(obj,mobj,idx);
            if isempty(bvobj), return; end
            %extract sampling box dimensions used for the volume calcs
            metatxt = bvobj.Data.Dataset.MetaData; 
            mxmn = getRangefromString(obj,metatxt);

            [K,dfv,ipr,bmv] = Sim_BMVfitting(bmv,prwvs,prdst,mxmn);
            if isempty(K)
                return; 
            elseif size(K,1)>1    %NOT TESTED
                prtime = datetime(getabstime(tset.pvts));
                labtxt = sprintf('BMV fits for profile %s using %s',...
                                    prdst.Description,wvdst.Description);
                [a,b,r2] = plotBMVfitting(obj,prtime(ipr),dfv,K,labtxt);
                    %------
                temp = array2table(K);
                props = {'K1','K2','K3','K4','K5','K6','K7','K8'};
                temp.Properties.VariableNames = props;
                temp.Properties.RowNames = cellstr(prtime(ipr));
                % save('temp.mat','temp')
                answer = questdlg('Save adjusted coefficients?','BMV model','Save','Quit','Quit');
                if strcmp(answer,'Quit'), return; end
                % [a,b] = adjustDissipationCoefficients(obj,a,b,params);
                fitcoeffs = struct('a',a,'b',b,'r2',r2);
                %the following ONLY works if there are no common fields 
                %(duplicate fields) in the 2 structures.
                params = getRunProperties(bmv);
                mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],...
                                            [fieldnames(x);fieldnames(y)]);
                fitcoeffs = mergestructs(fitcoeffs,params);           
                %SAVE in CT_BeachAnalysis.FitCoeffs
                mobj.Cases.CT_BeachAnalysis(bvrec).FitCoeffs = fitcoeffs;
            end
            getdialog('Run complete');
        end
%%
        function obj = runSimulation(mobj)    %NOT TESTED
            %function to run a simple 2D diffusion model
            obj = Sim_BMV;                      
            dsp = modelDSproperties(obj);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj);             
            %--------------------------------------------------------------------------
            % Model code 
            %--------------------------------------------------------------------------
            obj.UIsel = gobj.UIselection;
            obj.UIset = gobj.UIsettings;
            
            %get input data
            idx = [];         %subset of variables to be checked
            [bvobj,prdst,wvdst,bmv,prwvs,~] = getSimData(obj,mobj,idx);
            if isempty(bvobj), return; end
            %extract sampling box dimensions used for the volume calcs
            metatxt = bvobj.Data.Dataset.MetaData; 
            mxmn = getRangefromString(obj,metatxt);

            %call model function
            results = runBMVsimulation(bmv,mnvals,mxmn,prdst);  
            if isempty(results), return; end
            V = results{1}; X1 = results{6}; Z1 = results{7};
            labtxt = sprintf('%s using %s',prdst.Description,wvdst.Description);
            plotBMVmodel(obj,prtime,V,X1,Z1,labtxt);
            answer = questdlg('Save the results?','BMV model','Save','Quit','Quit');
            if strcmp(answer,'Quit'), return; end
            %--------------------------------------------------------------------------
            % Assign model output to a dstable using the defined dsproperties meta-data
            %--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(V,Z1,X1,'RowNames',prdst.RowNames,'DSproperties',dsp);

            %data to be displayed on Calcs Tab (may need to change tab list)
            output = volumesModel(bvobj,mobj);
            if ~isempty(output)                        
                src = findobj(mobj.mUI.Tabs.Children,'Tag','Volumes'); 
                headtxt = sprintf('Profile: %s\nFor: %s',...
                              output.source,output.metatxt);
                obj.TabOutput.table = output.table;
                obj.TabOutput.headtxt = headtxt;
                tablefigure(src,headtxt,output.table);
            end 
            %--------------------------------------------------------------------------
            % Save results
            %---------------------------------------------------------------
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = mtxt;
            %save results
            setDataSetRecord(obj,muicat,dst,'model');
            getdialog('Run complete');             
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private)
        function [bvobj,prdst,wvdst,bmv,prewvs,bvrec] = getSimData(obj,mobj,idx)
            %extract the data selected by the user
            %returns the class object of the beach volume data and the
            %dstables of the beach profile, wave data preceeding waves at
            %times of profile surveys (prewvs)
            prewvs=[]; prdst=[]; wvdst = []; bmv=[];
            muicat = mobj.Cases;
            hw = waitbar(0, 'Loading data. Please wait');
            
            %retrieve selected beach volume data
            bvcaserec = obj.UIsel(2).caserec;
            [bvobj,bvrec,catrec] = getCase(muicat,bvcaserec);
            %check that selection is correct
            if ~isa(bvobj,'CT_BeachAnalysis'), return; end

            %retrieve the associated beach profile dataset
            casedesc = catrec.CaseDescription;
            pid = split(casedesc,'-'); pid = pid{2};
            prrec = find(strcmp(muicat.Catalogue.CaseDescription,pid));
            prdst = getDataset(muicat,prrec,1);  %idset = 1
            waitbar(0.2)
            
            %retrieve an inshore wave data set
            wvcaserec = obj.UIsel(1).caserec;
            wvobj = getCase(muicat,wvcaserec);
            %check that selection is correct (rather than call getDataset)
            if ~isa(wvobj,'ctWaveModel'), return; end
            wvdst = getWaveModelDataset(wvobj,mobj,'Inwave_model',{'Tp'},...
                                                             wvcaserec);  
            %get the overlapping portion for the two data sets
            [wvdst,bvdst] = getoverlappingtimes(wvdst,bvobj.Data.Dataset,true);
            if prdst.RowRange{1}~=bvdst.RowRange{1} || prdst.RowRange{2}~=bvdst.RowRange{2}
                prdst = getsampleusingtime(prdst,bvdst.RowRange{1},bvdst.RowRange{2});
            end
            waitbar(0.4)

            %initialise the Sim_BMVmodel class
            if isempty(idx)
                idx = bvobj.FitCoeffs;
            end 
            BMVparams = getClassObj(mobj,'Inputs','Sim_BMVinput');
            bmv = setBMVmodelParams(BMVparams,mobj,idx);
            waitbar(0.6)
            
            %get the wave parameters for the interval before each survey
            %using the duration defined by wave averaging interval (days)
            %****IF prewvs NOT NEEDED FOR SIMULATION PROVIDE OPTION TO
            %aVOID BECAUSE IT TAKES A LONG TIME            
            waveint = bmv.WaveAvInt;
            prewvs = getpreceedingdata(wvdst,bvdst,waveint);            
            waitbar(1)

            %mean(Hs),std(Hs),mean(Tp)used in option 3 of ClosureOption
            %for Hallermeier zones
            bmv.meanWaves.mH = mean(wvdst.Hsi);  
            bmv.meanWaves.sH = std(wvdst.Hsi);  
            bmv.meanWaves.mT = mean(wvdst.Tp);
            close(hw);
        end  
%%
        function mxmn = getRangefromString(~,metatxt)
            %parse text string of max and min values and return numeric values
            [~,stid] = regexp(metatxt,'from ');
            edid = regexp(metatxt,' to');
            sep = regexp(metatxt,';')-1;
            mxmn.xmin = str2double(metatxt(stid(1):edid(1)));
            mxmn.xmax = str2double(metatxt(edid(1)+4:sep));
            mxmn.zmin = str2double(metatxt(stid(2):edid(2)));  
            mxmn.zmax = str2double(metatxt(edid(2)+4:end));
        end 
%%
        function [a,b] = adjustDissipationCoefficients(~,a,b,params)
            %regression sugggests coefficients are function of dim. fall velocity
            %extract resultant coefficients based on formulation used to make fit
            if params.BMVcoeffOption==1
        %         a(4) = ;
        %         a(5) = ;
        %         a(6) = ;
        %         a(7) = ;
        %         b(4) = ;
        %         b(5) = ;
        %         b(6) = ;
        %         b(7) = ;
            elseif params.BMVcoeffOption==2

            else
                warndlg('Incorrect value of BMVcoeffOption in CT_SimUI.adjustDissipationCoefficients')
            end
        end
%%
        function [a,b,r2] = plotBMVfitting(obj,tb1,dfv,K,labels)
            %generate plots of fir parameters for the BMVmodel
            hf = figure('Name','BMVfitting','Tag','PlotFig','Units','normalized');
            %move figure to top right
            hf.Position(1) = 1-hf.Position(3)-0.01;
            hf.Position(2) = 1-hf.Position(4)-0.12;

            ncoeff = size(K,2);
            ns = (ncoeff)/2;
            a = zeros(ncoeff,1); b = a; r2 = a;
            
%             subplot(ns,2,1)
%             plot(tb1,dfv)
%             ylabel('dfv')
%             xlabel('Time')             
            for k = 1:ncoeff
                subplot(ns,2,k)
                plot(dfv,K(:,k),'o','DisplayName','Fit estimates')
                ylabel(sprintf('K(%d)',k))
                [a(k),b(k),r2(k)] = addRegressionLine(obj,dfv,K(:,k),'Linear',3);
                legend('Location','best')
                if ncoeff==7 && (k==6 || k==7)
                    xlabel('dfv')
                elseif ncoeff==3 && (k==2 || k==3)
                    xlabel('dfv')
                end
            end   
            sgtitle(labels)
        end    
%%
        function [a,b,r2] = addRegressionLine(~,depdata,inddata,model,nint)
            [a,b,r2,x,y,txt] = regression_model(depdata,inddata,model,nint);
            hold on
            plot(x,y,'DisplayName',txt)
            hold off
        end
%%
        function plotBMVmodel(~,tb1,v,xc,zc,labels)
            %generate plot of outputs from predictive BMVmodel
            hf = figure('Name','BMVmodel','Tag','PlotFig','Units','normalized');
            %move figure to top right
            hf.Position(1) = 1-hf.Position(3)-0.01;
            hf.Position(2) = 1-hf.Position(4)-0.12;

            sgtitle(labels)
            subplot(2,1,1)
            plot(tb1,v)
            ylabel('Volume')
            subplot(2,1,2)
            plot(xc,zc,'--+')
            xlabel('x-centroid') 
            ylabel('z-centroid')            
        end
%%
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                      
                'Name',{'Volume','Elevation','Chainage'},...
                'Description',{'Volume','Elevation','Chainage'},...
                'Unit',{'m^3','mOD','m'},...
                'Label',{'Volume (m^3)','Elevation (mOD)','Chainage (m)'},...
                'QCflag',repmat({'model'},1,3));
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time'},...
                'Format',{'dd-MM-yyyy HH:mm:ss'});        
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});  
        end
    end      
end