classdef Sim_BMVinput < muiPropertyUI 
%
%-------class help---------------------------------------------------------
% NAME
%   Sim_BMVinput.m
% PURPOSE
%   Class to handle input data for BMV fitting and simulation model
% USAGE
%   obj = Sim_BMVinput.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   Sim_BMVfitting.m, Sim_BMVmodel.m, CT_Simulation.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Sediment sphericity (-)',...
                          'Beach porosity (-)',...
                          'Sediment grain size (m)',...                         
                          'Upper dune slope (1:uds)',...
                          'Upper beach slope (1:ubs)',...
                          'Beach crest, Zc or [Xc,Zc] (m,mOD)',...                          
                          'High water level (mOD)',...
                          'Low water level (mOD)',...
                          'Minimum distance to high water (m)',...
                          'Setup coefficient (-)',...
                          'Runup coefficient (default=1)',...
                          'Dune berm factor (0-1)',...
                          'Depth limited breaker index',...
                          'Wave averaging interval (days)',...
                          'Empirical coefficient option (1 or 2)',...
                          'Profile zones equation option (1,2,3)',...
                          'Beach adjustment method (1 or 2)', ...                          
                          'Fit using 3,4,6 or 8 parameters (3,4,6,8)'}
            %abstract properties in PropertyInterface for tab display
            TabDisplay   %structure defines how the property table is displayed
    end
    
    properties
        SedSphericity = 0.9;        %sphericity of sediment
        BeachPorosity = 0.35;       %beach porosity
        GrainSize = 0.001           %sediment grain size (m)       
        UpperDuneSlope              %upper dune slope (1:ds)
        UpperBeachSlope             %bed slope (1:bs)                    
        BeachCrest                  %beach crest, Zc or [Xc,Zc] (m,mOD)
        HWlevel                     %high water level (mOD)
        LWlevel                     %low water level (mOD)
        xHWmn                       %minimum distance to high water (m)
        SetupCoeff = 0.25;          %setup coefficient
        RunupCoeff = 1.0;           %runup coefficient
        DuneBermFactor = 0.3;       %dune berm factor  
        BreakerIndex = 0.78;        %Depth limited breaker index
        WaveAvInt = 2;              %Wave averaging interval (days)
        BMVcoeffOption = 1;         %Empirical coefficient option (1 or 2)
        ClosureOption = 1;          %Profile zones equation option (1,2,3)
        AdjustMethod = 1;           %Beach adjustment method (1 or 2)        
        BMVfitOption = 6;           %Use 3,4,6 or 8 fit parameters (1 or 2)
    end    
  
%%   
    methods (Access=protected)
        function obj = Sim_BMVinput(mobj)  
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI fcn
            
            %to use non-numeric entries then one can either pre-assign 
            %the values in properties defintion, above, or specity the 
            %PropertyType as a cell array, e.g.:
            % obj.PropertyType = [{'datetime','string','logical'},...
            %                                       repmat({'double'},1,8)];
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'Sim_BMVinput';  
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = Sim_BMVinput(mobj);            
            end
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end     
    end
%%        
    %add other functions to operate on properties as required
    methods
        function pfit = selectFitParams(obj,fitcoeffs)
            %allow user to edit the fit coefficients for use in BMV model
            %fitcoeffs is a struct of regression coefficients and R2
            %values created in Sim_BMVfitting and set in CT_Simulation
            if isempty(fitcoeffs)
                %edit static input parameters based on subset defined by idx
                idx = 12:14;         %subset of variables to be checked
                editPropertySubset(obj,idx);
                params = getPropertiesStruct(obj);
                return;
            else
                params = getPropertiesStruct(obj);
            end
            f = fitcoeffs;
            %allow user to adjust the input values
            nprop = 7;  allprop = nprop+2;
            defvals = cell(allprop,1);
            prompt = {'K1 (a+b.x)','K2 (a+b.x)','K3 (a+b.x)','K4 (a+b.x)',...
                'K5 (a+b.x)','K6 (a+b.x)','K7 (a+b.x)',...
                'Beach adjustment method (1 or 2)', ...
                'Empirical coefficient option (0, 1 or 2)'};
            for i=1:nprop
                defvals{i} = sprintf('%.4f   %.4f',f.a(i),f.b(i));
                prompt{i} = sprintf('%s: R2 = %.3f',prompt{i},f.r2(i));
            end
            defvals{8} = num2str(f.AdjustMethod);
%             defvals{9} = num2str(f.BMVcoeffOption);
            defvals{9} = '0';
            answer = inputdlg(prompt,'BMV coefficients',1,defvals);
            if isempty(answer), return; end
            
            %assign the selcted values to the params struct
            props = {'Kc1','Kc2','Kc3','Kc4','Kc5','Kc6','Kc7'};
            for i=1:nprop
                pfit.(props{i}) = str2num(answer{i}); %#ok<ST2NM> - handles arrays str2double does not
            end
            
%             params.AdjustMethod = str2double(answer{8});
%             params.BMVcoeffOption = str2double(answer{9});
%             params.WaveAvInt = f.WaveAvInt;
%             
%             %the following ONLY works if there are no common fields
%             %(duplicate fields) in the 2 structures.
%             mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],...
%                 [fieldnames(x);fieldnames(y)]);
%             params = mergestructs(params,pfit);
            
        end
%%
        function bmv = setBMVmodelParams(obj,mobj,idx)
            %assign Sim_BMVinput properties and instantiate an instance
            % bmv - BMVmodel class
            % obj - instance of Sim_BMVinput
            % idx - if fitting idx is subset of parameters to adjust
            %       if simulation idxis struct of fit coefficients
            %instantiate BMVmodel
            bmv = Sim_BMVmodel();

            %get BMVmodel fit coefficients
            if isstruct(idx)
                %edit dynamic fit parameters and add fitting coefficients
                %assigned to input volume property FitCoeffs to 'params'
                % - additional coeffs allows functions of dfv to be used in model
                bmv.Coeffs = selectFitParams(obj,prv.obj.FitCoeffs);
            else
                %edit static input parameters based on subset defined by idx
                editPropertySubset(obj,idx);
                bmv.Coeffs = [obj.SetupCoeff,obj.RunupCoeff,...
                                    obj.DuneBermFactor,obj.BreakerIndex];
            end

            bmvprops = {'S','phi','d50','ds','ubs','XZc','Zhw','Zlw','xHWmn'};
            bmvinput = {'SedSphericity','BeachPorosity','GrainSize',...
                        'UpperDuneSlope','UpperBeachSlope','BeachCrest',...
                        'HWlevel','LWlevel','xHWmn'};
            %transfer input properties to model        
            for k=1:length(bmvprops)
                bmv.(bmvprops{k}) = obj.(bmvinput{k});
            end

            %add system variables
            bmv.g = mobj.Constants.Gravity;               %gravity (m/s^2)
            bmv.visc = mobj.Constants.KinematicViscosity; %viscosity (m2/s)
            bmv.rhow = mobj.Constants.WaterDensity;       %water density (kg/m3)
            bmv.rhos = mobj.Constants.SedimentDensity;    %sediment density (kg/m3)

            %set the run time properties  
            vals = [obj.WaveAvInt,obj.BMVcoeffOption,obj.ClosureOption,...
                                        obj.AdjustMethod,obj.BMVfitOption];
            bmv = setRunProps(bmv,vals);
        end
    end
end