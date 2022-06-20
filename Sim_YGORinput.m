classdef Sim_YGORinput < muiPropertyUI 
%
%-------class help---------------------------------------------------------
% NAME
%   Sim_YGORinput.m
% PURPOSE
%   Class to handle input data for YGOR fitting and simulation model
% USAGE
%   obj = Sim_YGORinput.setParamInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   Sim_YGORmodel.m, CT_Simulation.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2019
%--------------------------------------------------------------------------
%
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Forecast origin position (m):',...
                          'Forecast origin date (dd-MM-yyyy):',...
                          'Equilibirum slope, a:',...
                          'Equilibirum intercept, b:',...
                          'Accretion rate coefficeint, C+:',...
                          'Erosion rate coefficient, C-:',...
                          'Offset to initial shore position:',...
                          'Fitted coefficients RMSE (read only)'};
        %abstract properties in PropertyInterface for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        ForecastOriginPosition = 0        %Forecast origin position (m):
        ForecastOriginDate ='01-01-1990'  %Forecast origin date (dd-mmm-yyyy)
        EqSlope = 0                       %Equilibirum slope
        EqIntercept = 0                   %Equilibirum intercept
        AccretionRateCoeff = 0            %Accretion rate coefficeint
        ErosionRateCoeff = 0              %Erosion rate coefficient
        InitialOffset = 0                 %Offset to initial shore position
        YGORrmse = 0                      %RMSE for fitted coefficients (read only)
    end    

%%   
    methods (Access=?Sim_YGOR)
        function obj = Sim_YGORinput(mobj)   
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
            classname = 'Sim_YGORinput';      
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = Sim_YGORinput(mobj);            
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
        function simparams = getSimParams(obj)
            %get the structure properties as a struct
            sp.Sh0  = obj.ForecastOriginPosition; %Forecast origin position (m)
            sp.Dat0 = obj.ForecastOriginDate;     %Forecast origin date (dd-mmm-yyyy)
            sp.C(1) = obj.EqSlope;                %Equilibirum slope
            sp.C(2) = obj.EqIntercept;            %Equilibirum intercept
            sp.C(3) = obj.AccretionRateCoeff;     %Accretion rate coefficient
            sp.C(4) = obj.ErosionRateCoeff;       %Erosion rate coefficient
            sp.C(5) = obj.InitialOffset;          %Offset to initial shore position
            simparams = sp;
        end
    end          
end