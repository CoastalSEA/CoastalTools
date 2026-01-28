classdef CT_WaveModels < muiDataSet                  
%
%-------class help------------------------------------------------------
% NAME
%   Model_template.m
% PURPOSE
%   Class for Model CT_WaveModels to be run from CoastalTools
%
% SEE ALSO
%   muiDataSet
%   For full list of functions used:
%   flist = matlab.codetools.requiredFilesAndProducts('CT_BeachAnalysis.m')'
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties: 
    end
    
    properties (Hidden)
        ModelType        %model used for the particular instance
    end
    
    properties (Transient)
        MenuList = {'Littoral Drift','X-shore Transport','Wave Power',...
                     'Runup','Structure Overtopping','Beach Overtopping',...
                     'Iribarren Number','Beach type'}
        ModelName = {'Drift','Xshore','WavePower','Runup','StructureOtop',...
                                      'BeachOtop','Iribarren','BeachType'};
    end
    
    methods
        function obj = CT_WaveModels()                    
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj,src)
            %function to run a simple 2D diffusion model
            obj = CT_WaveModels;    
            id_model = find(strcmp(obj.MenuList,src.Text));            
            dsp = modelDSproperties(obj,id_model);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            %sturct input allows class to be used as well as ctWaveModel
            txt1 = 'Use Setup to define model input parameters';
            txt2 = 'Run Nearshore wave model to define input waves';
            if id_model==4
                altclass = struct('idv',{2},'class',{'ctWaveData'});
                isvalid = isValidModel(mobj, metaclass(obj).Name,altclass);
            else
                isvalid = isValidModel(mobj, metaclass(obj).Name);
            end
            if ~isvalid, warndlg(sprintf('%s\n%s',txt1,txt2)); return; end              

            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); %can be updated in model calls (eg drift)
%--------------------------------------------------------------------------
% Model code
%--------------------------------------------------------------------------
            site = mobj.Inputs.ctWaveParameters;  
            switch id_model
                case 1              %Littoral Drift
                    output = driftModel(obj,mobj,site);
                case 2              %Cross-shore transport
                    output = xshoreModel(obj,mobj,site);
                case 3              %Wave Power
                    output = energyModel(obj,mobj,site);
                case 4              %Runup                    
                    output = runupModel(obj,mobj,site);
                case 5              %Strucure Overtopping
                    output = sOtopModel(obj,mobj,site);
                case 6              %Strucure Overtopping
                    output = bOtopModel(obj,mobj,site);
                case 7              %Iribarren Number
                    output = iribarrenModel(obj,mobj,site); 
                case 8              %Beach type based on fall velocity
                    output = beachTypeModel(obj,mobj,site);
            end
            if isempty(output) || ~isfield(output,'results') || isempty(output.results{1})
                return; %user cancelled or no results returned
            end 
            %setRunParam(obj,mobj,output.wvrec); 
            setRunParam_InputDataSet(obj,mobj,output.wvrec)
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(output.results{:},'RowNames',output.modeltime,...
                                                    'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            obj.ModelType = obj.ModelName{id_model};
            dst.Source = sprintf('Class %s, using %s model',metaclass(obj).Name,...
                                                            obj.ModelType);
            dst.MetaData = output.metatxt;
            %save results
            setDataSetRecord(obj,muicat,dst,'model');
            getdialog('Run complete');
        end
    end
%%
    methods
        function output = driftModel(obj,mobj,site)
            % compute the littoral drift based on wave heights at the edge
            % of the surf zone
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by model:
                % Hsi  - inshore significant wave height (m)
                % Tp   - peak wave period (s)
                % Diri - wave direction (degrees TN)
                % depi  - water depth at wave point (m)
                % theta- angle of shoreline from north (degrees TN)                
                % d50  - grain size d50 (m)
                % Kc   - drift coefficient (-)
            % FUNCTION CALLS (external) - littoraldrift

            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = []; return; end
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data
            
            %bed slope within surf zone (half depth of inshore wave point)  
            ubs = site.UpperBeachSlope;
            z1km = site.BedLevelat1km;            
            bs = profileslope(wv.depi/2,wv.swl,z1km,ubs); %first argument is depth
            
            %call drift model and add longhshore drift to wave time series
            g = mobj.Constants.Gravity;
            rhw = mobj.Constants.WaterDensity;
            rhs = mobj.Constants.SedimentDensity;
            vsc = mobj.Constants.KinematicViscosity;
            [theta,d50,Kc] = getDriftSettings(obj,site);

            Qall = littoraldrift(wv.Hsi,wv.Tp,wv.Diri,wv.depi,...
                                        theta,bs,d50,Kc,g,rhs,rhw,vsc);            
            %user selects which model to use
            tlist = {'CERC formula (SPM, 1994)',...
                     'Dynamics of Marine Sands, Soulsby',...
                     'Kamphuis formula',...
                     'Damgaard & Soulsby (shingle)'};
            [h_dlg,ok] = listdlg('Name','Plot profile', ...
                                 'PromptString','Select formula', ...
                                 'ListSize',[200,80], ...
                                 'SelectionMode','single', ...
                                 'ListString',tlist);
            if ok==0, return; end 
            Q = Qall(:,h_dlg);
            formula = tlist{h_dlg};
            bs = mean(bs,'omitnan');
            zi = mean((wv.swl-wv.depi),'omitnan');
            
            output.results = {Q};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('Drift using %s; Theta=%g; d50=%g; Kc=%g; Beach slope=1:%.1f; Zi=%g',...
                                    formula,theta,d50,Kc,bs,zi);    
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);    
        end

%%
        function output = xshoreModel(~,mobj,site)
            % compute the xross-shore transport based in the surf zone
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by model:
                % Hsi  - inshore significant wave height (m)
                % Tp   - peak wave period (s)
                % Diri - wave direction (degrees TN)
                % depi  - water depth at wave point (m)
                % theta- angle of shoreline fron north (degrees TN)                
                % d50  - grain size d50 (m)
                % Kc   - drift coefficient (-)
            % FUNCTION CALLS (external)
            % xshore_bailard
            
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = []; return; end           
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data
            
            %bed slope within surf zone (half depth of inshore wave point)  
            ubs = site.UpperBeachSlope;
            z1km = site.BedLevelat1km;            
            bs = profileslope(wv.depi/2,wv.swl,z1km,ubs); %first argument is depth
            
            %call xshore model and add longhshore drift to wave time series
            g = mobj.Constants.Gravity;
            rhw = mobj.Constants.WaterDensity;
            rhs = mobj.Constants.SedimentDensity;
            vsc = mobj.Constants.KinematicViscosity;
            theta = site.ShorelineAngle;
            d50 = site.GrainSize;
            Qx = xshore_bailard(wv.Hsi,wv.Tp,wv.Diri,wv.depi,...
                                            theta,bs,d50,g,rhw,rhs,vsc);
            
            bs = mean(bs,'omitnan');
            zi = mean((wv.swl-wv.depi),'omitnan');
            
            output.results = {Qx};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('d50=%g; Beach slope=1:%.1f; Zi=%g',...
                                           d50,bs,zi);
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
        end

%%
        function output = energyModel(~,mobj,site)
            %calculate the inshore wave power (wave energy flux) using 
            %linear wave theory and the refracted, shoaled and 
            %breaking wave conditions at the defined inshore point
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by model:
                % swl - still water level (m above datum)
                % Hs  - incident significant wave height (m)
                % Tp  - peak wave period (s)
                % zi   - inshore bed level (mOD) NaN for edge of surf zone
                % ubs  - upper beach slope (1:ms)
                % z1km - bed elevation 1km from shoreline (mOD)
                % bl  - bed level(m above datum)
                % rhow- water density(kg/m3)
            % FUNCTION CALLS (external)
            % wave_energyflux
            
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = [];  return; end
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data
            
            %average nearshore bed slope - for meta data only
            zi = mean((wv.swl-wv.depi),'omitnan');
            dep = mean((wv.depi),'omitnan');
            ubs = site.UpperBeachSlope;
            z1km = site.BedLevelat1km;            
            bs = profileslope(dep/2,wv.swl,z1km,ubs);
            
            %calculate energy flux and add results to wave timeseries
            rhow = mobj.Constants.WaterDensity;
            g = mobj.Constants.Gravity;
            Hrms = wv.Hsi/sqrt(2);            %rms wave height
            Ef = wave_energyflux(Hrms,wv.Tp,wv.depi,rhow,g);
 
            output.results = {Ef};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('Nearshore slope=1:%.1f; Zi=%g',bs,zi);   
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
        end

%%
        function output = runupModel(obj,mobj,site)
             %calculate the 2% runup magnitude and elevation (R2+WL) using
            %the equation of Stockdon etal, 2006. Uses the effective
            %offshore wave conditions by shoaling the measured data to deep
            %water and the beach foreshore slope (slope at water line).
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by model:
                % swl - still water level (m above datum)
                % Hs  - measured (offshore) significant wave height (m)
                % Tp  - peak wave period (s)
                % dep0 - measurement site (offshore) water depth (m)
                % ubs  - upper beach slope (1:ms)
                % z1km - bed elevation 1km from shoreline (mOD)
                % visc - viscosity (m/s^2)
                % d50  - sediment median grain size (m)
                % S    - sediment sphericity (-)
                % phi  - porosity of sand bed (-) 
            % FUNCTION CALLS (external)
            % shoaling, runup, iribarren, runup_slope          
            
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            [inwave,~,~] = selectCaseObj(mobj.Cases,[],...
                    {'ctWaveModel','WRM_WaveModel','ctWaveData'},promptxt);
            if isempty(inwave), output = [];   return; end            
            %retrieve an inshore wave data set
            if isa(inwave,'ctWaveData')
                %retrieve an offshore (e.g. buoy) wave data set
                % %not implemented
                % wv = [];
                % msgbox('Select nearshore (model) wave data')
                output.wvrec = caseRec(mobj.Cases,inwave.CaseIndex);
                [wv,meta] = addWaveWLdataset(inwave,mobj,output.wvrec);
                mtxt2 = meta.inptxt;
                idv = strcmp(wv.VariableNames,'Hs');
                wv.VariableNames{idv} = 'Hsi';%match variable names to wave model
            else
                %retrieve an inshore wave data set
                [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
                mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            end
            if isempty(wv), output = []; return; end %user cancelled or no data
            
            %get "effective" deepwater offshore wave heights   
            if any(strcmp(wv.VariableNames,'depi'))
                dep0 = wv.depi;        %depth at inshore point used by wave model
            else
                dep0 = wv.swl-site.OffshoreBedLevel;
                wv = addvars(wv,dep0,'NewVariableNames','depi');
            end
            dep1  = 100;           %offshore deep water depth           
            Hs0 = shoaling(wv.Hsi,wv.Tp,dep0,dep1);
            
            %find slope at swl on beach  
            swl_0  = wv.swl; 
            swl_0(isnan(swl_0)) = 0; %replace nans with zero
            ubs = site.UpperBeachSlope;
            z1km = site.BedLevelat1km;            
            bs = profileslope(0,swl_0,z1km,ubs); %first argument is depth
            
            %get runup and elevation of runup             
            [R2,~,~,gR2] = runup(bs,Hs0,wv.Tp);            
            zR2 = R2+wv.swl;
            zgR2 = gR2+wv.swl;
            
            %get runup slope using method of Reis and Gama, 2010
            Rslope = getRunupSlope(obj,mobj,Hs0,wv.Tp,site);
            
            bs = mean(bs,'omitnan'); %use average slope in metadata
            zi = mean((wv.swl-wv.depi),'omitnan');
            
            output.results = {R2,zR2,gR2,zgR2,Rslope};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('Beach slope=1:%0.1f; mean zi=%g',bs,zi); 
            
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
        end

%%
        function output = sOtopModel(~,mobj,site)
            %calculate the wave overtopping for the defined structure
            %properties using the HRW Owen model.
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by overtopping model:
                % Hsi  - inshore significant wave height (m)
                % Tp   - peak wave period (s)
                % Diri - wave direction (degrees TN)
                % swl  - still water level (m above datum)
                % alp  - wave angle in degrees from normal to the wall
                %      - taken to be the same as the shoreline angle
                % zi   - inshore bed level (mOD) NaN for edge of surf zone
            % structure definition (obj.Structure)
            % FUNCTION CALLS (external)
            % otop_Q (uses celerity, hb_break, coeff_AB)

            if isfield(mobj.Inputs,'ctStructureInput') && ...
                                    ~isempty(mobj.Inputs.ctStructureInput)
                structprops = mobj.Inputs.ctStructureInput;   
            else
                warndlg('Structure properties have not been defined')
                output = [];
                return;
            end
            
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = [];  return; end 
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data
                                       
            varnames = wv.VariableNames;
            if ~any(strcmp(varnames,'Tz'))
                Tz = wv.Tp*0.7775;  %assumes JONSWAP with gamma=3.3
                wv = addvars(wv,Tz,'NewVariableNames','Tz'); 
                getdialog('Tz estimated using Tp')
            end
            
            %bed slope half wave length in front of structure
            ubs = site.UpperBeachSlope;
            z1km = site.BedLevelat1km;       
            bs = profileslope(wv.depi/2,wv.swl,z1km,ubs); %first argument is depth

            g = mobj.Constants.Gravity;
            theta = site.ShorelineAngle;
            alp = wv.Diri-(theta+90);            
            structure = getStructure(structprops);
            Q = otop_Q(wv.swl,wv.Hsi,wv.Tz,alp,bs,g,structure);            
            bs = mean(bs,'omitnan'); %use average slope in metadata

            output.results = {Q};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('Neashore slope=1:%0.1f; theta=%g',bs,theta); 
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
        end

%%
        function output = bOtopModel(~,mobj,site)
            %calculate the wave overtopping for a gravel beach
            %properties required by the model of Stoke's et al, 2021
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % site - uses ctWaveParameters properties
            %        BeachCrestLevel  - crest elevation (mOD)
            %        UpperBeachSlope - upper beach slope (1:ubs)
            % Variables used by overtopping model:
                % Hsi  - inshore significant wave height (m)
                % Tp   - peak wave period (s)
                % swl  - still water level (m above datum)
            % FUNCTION CALLS (external)
            % otopBeach (uses hb_break, profileslope)
            
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = [];  return; end 
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data

            %get "effective" deepwater offshore wave heights  
            if any(strcmp(wv.VariableNames,'depi'))
                dep0 = wv.depi;        %depth at inshore point used by wave model
            else
                dep0 = wv.swl-site.OffshoreBedLevel;
            end
            dep1  = 100;           %offshore deep water depth           
            H0 = shoaling(wv.Hsi,wv.Tp,dep0,dep1);

            inp = inputdlg('Beach toe level (mOD):','Beach toe',1,{'0'});
            if isempty(inp), output = []; return; end    
            beach = getPropertiesStruct(site);
            beach.BeachToeLevel = str2double(inp{1});

            [Q,gR2,zgR2] = otopBeach(wv.swl,wv.Hsi,wv.Tp,H0,beach);

            output.results = {Q,gR2,zgR2};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('Beach toe level %g',beach.BeachToeLevel); 
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
        end

%%
        function output = iribarrenModel(~,mobj,site)
            %calculate the Iribarren number and breaker type at the 
            %inshore point using incident wave conditions (not Ho).
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by model:
                % swl - still water level (m above datum)
                % Hs  - incident significant wave height (m)
                % Tp  - peak wave period (s)
                % zi  - inshore bed level (mOD) NaN for edge of surf zone
                % ubs  - upper beach slope (1:ms)
                % z1km - bed elevation 1km from shoreline (mOD)
            % FUNCTION CALLS (external)
            % Iribarren (uses Hs_break, Hb_break)
            
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:'; 
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = [];  return; end         
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data
            
            %find slope at swl on beach         
            ubs = site.UpperBeachSlope;
            z1km = site.BedLevelat1km;            
            bs = profileslope(0,wv.swl,z1km,ubs); %first argument is depth
            
            %call Iribarren function and add results to wave timeseries
            g = mobj.Constants.Gravity;
            [Iri,Ityp] = iribarren(wv.Hsi,wv.Tp,bs,g);   
            bs = mean(bs,'omitnan');
            zi = mean((wv.swl-wv.depi),'omitnan');
            
            output.results = {Iri,Ityp};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('Beach slope=1:%.1f; Zi=%g',bs,zi); 
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
        end

%%
       function output = beachTypeModel(~,mobj,site)
            % compute the beach type based on dimensionless fall velocity
            % INPUTS
            % mobj - handle to CoastlTools class objects
            % Variables used by model:
                % Hsi - incident significant wave height (m)
                % Tp  - peak wave period (s)
                % rhow- water density(kg/m3)
                % rhos- sediment density(kg/m3)   
                % visc- kinematic viscosity (m2/s)
                % d50 - sediment grain size(m)
                
            %class instance for inshore wave data
            promptxt = 'Select nearshore wave data set:';             
            inwave = selectCaseObj(mobj.Cases,[],{'ctWaveModel','WRM_WaveModel'},promptxt);
            if isempty(inwave), output = [];  return; end
            %retrieve an inshore wave data set
            [wv,output.wvrec] = getWaveModelDataset(inwave,mobj,...
                                                {'Inwave_model'},{'Tp'});
            if isempty(wv), return; end %user cancelled or no data

            %calculate dimensionless fall velocity and add results to wave timeseries
            d50 = site.GrainSize;
            g = mobj.Constants.Gravity;
            rhow = mobj.Constants.WaterDensity;
            rhos = mobj.Constants.SedimentDensity;
            visc = mobj.Constants.KinematicViscosity;
            ws = settling_velocity(d50,g,rhow,rhos,visc);  %sediment fall velocity
            dfv = wv.Hsi./ws./wv.Tp;  %dimensionless fall velocity
            %if dfv(i)>=6  'dissipative';
            %elseif dfv(i)<6 && dfv(ij)>1  'intermediate';
            %else  'reflective';
            output.results = {dfv};
            output.modeltime = wv.RowNames;
            mtxt1 = sprintf('d50=%g',d50); 
            mtxt2 = sprintf('Using %s case for wave input',wv.Description);
            output.metatxt = sprintf('%s\n%s',mtxt1,mtxt2);  
       end

%%
        function Rslope = getRunupSlope(~,mobj,Hs0,Tp,site)
            %setup input and call runup_slope
            g    = mobj.Constants.Gravity;            %gravity (m/s^2)
            visc = mobj.Constants.KinematicViscosity; %viscosity (m2/s)
            ubs  = site.UpperBeachSlope; %typical/average beach slope (1:bs)
            d50  = site.GrainSize;       %sediment grain size (m)         
            S    = 0.9;    %sediment sphericity, using mean sphericity of objects from cube to icosahedron
            P    = 0.35;   %porosity of sand bed (Soulsby, 1991)
            C    = 0.5;    %Setup coefficient (scales setup to estimate flow depth)
            prompt = {'Setup-coefficient','Sediment sphericity','Bed porosity'};
            defvals = {'0.5','0.9','0.35'};
            answer = inputdlg(prompt,'Runup slope',1,defvals);
            if ~isempty(answer)
                C = str2double(answer{1});
                S = str2double(answer{2});
                P = str2double(answer{3});
            end
            Rslope = runup_slope(Hs0,Tp,ubs,C,d50,S,P,g,visc);                        
        end 

%%
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private) 
        function [theta,d50,Kc] = getDriftSettings(obj,site)
            %check that the current sit parameters settings are correct
            %modifications used to update RunParams stored with Case
            theta = site.ShorelineAngle;
            d50 = site.GrainSize;
            Kc = site.DriftCoefficient;
            promptxt = {'Shoreline angle (degTN)','Sediment grain size (m)',...
                        'CERC Drift coefficient, Kc'};
            defaults = {num2str(theta),num2str(d50),num2str(Kc)};
            data = inputdlg(promptxt,'Drift settings',1,defaults);
            if isempty(data), return; end  %no change to default settings
            theta = str2double(data{1});
            d50 = str2double(data{2});
            Kc = str2double(data{3}); 
            %update RunParam settings
            obj.RunParam.ctWaveParameters.ShorelineAngle = theta;
            obj.RunParam.ctWaveParameters.GrainSize = d50;
            obj.RunParam.ctWaveParameters.DriftCoefficient = Kc;
        end

%%
        function dsp = modelDSproperties(~,id_model) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            switch id_model
                case 1             %Drift
                    dsp.Variables = struct(...                       
                        'Name',{'Qs'},...
                        'Description',{'Alongshore drift rate potential'},...
                        'Unit',{'m^3/s'},...
                        'Label',{'Transport rate (m^3/s)'},...
                        'QCflag',{'model'});
                case 2             %X-shore
                    dsp.Variables = struct(...                       
                        'Name',{'Qx'},...
                        'Description',{'Cross-shore transport rate'},...
                        'Unit',{'m^3/s'},...
                        'Label',{'Transport rate (m^3/s)'},...
                        'QCflag',{'model'});
                case 3             %Energy flux/wave power
                    dsp.Variables = struct(...                       
                        'Name',{'Eflux'},...
                        'Description',{'Inshore wave power'},...
                        'Unit',{'J/ms'},...
                        'Label',{'Wave Power (J/ms)'},...
                        'QCflag',{'model'});
                case 4              %Runup
                    dsp.Variables = struct(...                       
                        'Name',{'R2','zR2','gR2','zgR2','Rslope'},...
                        'Description',{'Runup (R2%)','Runup elevation (R2%+swl)',...
                                       'Runup gravel (gR2%)','Runup gravel elevation (gR2%+swl)',...
                                       'Runup slope (1:rs)'},...
                        'Unit',{'m','mOD','m','mOD','-'},...
                        'Label',{'Runup distance (m)','Runup elevation (mOD)',...
                                 'Runup distance (m)','Runup elevation (mOD)',...
                                 'Runup slope (1:rs)'},...
                        'QCflag',repmat({'model'},1,5));
                case 5              %Structure Overtopping
                    dsp.Variables = struct(...                       
                        'Name',{'Qotop'},...
                        'Description',{'Overtopping discharge'},...
                        'Unit',{'m^3/s'},...
                        'Label',{'Overtopping discharge (m^3/s)'},...
                        'QCflag',{'model'});
                case 6
                    dsp.Variables = struct(...                       
                        'Name',{'Qotop','gR2','zgR2'},...
                        'Description',{'Overtopping discharge',...
                                       'Runup gravel (gR2%)',...
                                       'Runup gravel elevation (gR2%+swl)'},...
                        'Unit',{'m^3/s','m','mOD'},...
                        'Label',{'Overtopping discharge (m^3/s)',...
                                 'Runup distance (m)','Runup elevation (mOD)'},...
                        'QCflag',repmat({'model'},1,3));
                case 7              %Iribarren
                    dsp.Variables = struct(...                       
                        'Name',{'Iri','BreakerType'},...
                        'Description',{'Iribarren wave number','Breaker type'},...
                        'Unit',{'-','-'},...
                        'Label',{'Iribarren wave number','Breaker Type (1-Spill; 2-Plunge; 3-Surge)'},...
                        'QCflag',{'model','model'});
                case 8              %Beach tupe 
                    dsp.Variables = struct(...                       
                        'Name',{'dfv'},...
                        'Description',{'Dimensionless fall velocity'},...
                        'Unit',{'-'},...
                        'Label',{'Dimensionless fall velocity'},...
                        'QCflag',{'model'});
            end
            %
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