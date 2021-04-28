classdef Sim_BMVmodel < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   Sim_BMVmodel.m
% PURPOSE
%   Class for generataion of a BMV equilibirum beach profile and
%   - fitting profiles to surveyed profile data
%   - simulating the profile change using a wave timeseries  
% SEE ALSO
%   Sim_BMVfitting.m, Sim_BMVinput.m, CT_Simulation.m
%   flist = matlab.codetools.requiredFilesAndProducts('Sim_BMVmodel.m')'
%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%    
    properties (Hidden)
        %abstract properties for a PropertyInterface that are not used here
        %definition of labels to be used for tabular display
        PropertyLabels = {'Sediment sphericity (-)',...
                          'Beach porosity (-)',...
                          'Sediment grain size (m)',...                          
                          'Upper dune slope (1:uds)',...
                          'Upper beach slope (1:ubs)',...
                          'Beach crest level (mOD)',...
                          'High water level (mOD)',...
                          'Low water level (mOD)',...  
                          'Minimum distance to high water (m)',...
                          'Acceleration due to gravity (m/s2)',...
                          'Water density (kg/m3)',...
                          'Sediment density (kg/m3)',...
                          'Kinematic viscosity (m2/s)'}
        TabDisplay       %structure defines how the property table is displayed
     end
    
    properties (Transient)
        BMVrunLabels = {'Wave averaging interval (days)',...
                        'Empirical coefficient option (1 or 2)',... 
                        'Closure depth equation option (1 or 2)',...
                        'Beach adjustment method (1 or 2)', ...                                             
                        'Fit using 3,4,6 or 8 parameters (3,4,6,8)'}
        BMVrunProps = {'WaveAvInt','BMVcoeffOption','ClosureOption'...
                                        'AdjustMethod','BMVfitOption'}          
        WaveAvInt = 2;              %Wave averaging interval (days)
        BMVcoeffOption = 1;         %Empirical coefficient option (1 or 2)
        ClosureOption = 1;          %Profile zones equation option (1 or 2)
        AdjustMethod = 1;           %Beach adjustment method (1 or 2)        
        BMVfitOption = 7;           %Use 3, 5 or 7 fit parameters (1 or 2)
        Coeffs                      %fit coefficient vector
        ismodel                     %true if simulation is being run
        Xref                        %X co-ordinate for reference profile
        Zref                        %Z co-ordinate for reference profile                    
        Zc                          %Z position of dune crest
        meanWaves                   %array of mean and std for Hs and Tp
    end

    properties
        S               %sphericity of sediment
        phi             %beach porosity 
        d50             %sediment grain size (m)
        ds              %upper dune slope (1:ds)
        ubs             %typical/average beach slope (1:bs)
        XZc             %crest position (mOD)      
        Zhw             %high water level (mOD)
        Zlw             %low water level (mOD)
        xHWmn           %minimum distance to high water over time could make this vectore for (a+b*t)
%         BMVcoeffOption  %Empirical coefficient option (1 or 2)
        %add system variables
        g = 9.81        %gravity (m/s^2)
        rhow = 1025     %water density (kg/m3)
        rhos = 2650     %sediment density (kg/m3)
        visc = 1.36e-6  %viscosity (m2/s)       
    end
% %%    
%     methods (Access=private)
%         function obj = Sim_BMVmodel()
%             %constructor to initialise object
%         end
%     end
% %%
%     methods (Static)
    methods
        function obj = Sim_BMVmodel(params)
            %constructor to initialise object
            %if params empty calls gui for user to set Property values
            %             obj = Sim_BMVmodel();
            if nargin<1
                return; %return an empty instance
            elseif isempty(params)
                nrec = 9;
                obj = editProperties(obj,nrec);  %add nrec to limit length of UI
            else
                obj = setProperties(obj,struct2cell(params));
            end
        end
%     end
% %%
%     methods
        function propdata = getBMVmodel(obj)
            %returns properties of BMVmodel as a struct
            propdata = getPropertiesStruct(obj);
        end     
%%
        function obj = setBMVmodel(obj,params,idx)
            %updates properties of BMVmodel with the values in params struct
            %if params is empty then uses the inpudlg to input values
            if nargin<2 || isempty(params)
                nrec = 9;
                obj = editProperties(obj,nrec);  %add nrec to limit length of UI    
            else
                if nargin<3
                    obj = setProperties(obj,struct2cell(params));
                else
                    obj = editPropertySubset(obj,idx);
                end
            end
        end
%%
        function obj = setRunProps(obj,vals)
            %set the run time properties
            %update with vals, if vals is empty prompt for values
            propnames = obj.BMVrunProps;
            nprop = length(propnames);
            if nargin<2 || isempty(vals)
                proptxt = obj.BMVrunLabels;
                defvals = cell(nprop,1);
                for k=1:nprop
                    defvals{k} = num2str(obj.(propnames{k}));
                end
                vals = inputdlg(proptxt,'Run time properties',1,defvals);
                if isempty(vals), return; end  %user cancelled                
                for k=1:nprop
                    obj.(propnames{k}) = str2double(vals{k});
                end                
            else
                for k=1:nprop
                    obj.(propnames{k}) = vals(k);
                end
            end
        end    
%%
        function obj = setCoeffProps(obj,vals)
            %set the coefficients used to define the first 3 values of
            %fitting vector K
            %update K with vals, if vals is empty prompt for values
            if nargin<2 || isempty(vals)
                proptxt = {'Setup coefficient (-)',...
                           'Runup coefficient (default=1)',...
                           'Dune berm factor (0-1)'};
                defvals = cell(3,1);
                for i=1:3
                    defvals{i} = num2str(obj.Coeffs(i));
                end
                vals = inputdlg(proptxt,'Model coefficients',1,defvals);
                if isempty(vals), return; end  %user cancelled
                for i=1:3
                    obj.Coeffs(i) = str2double(vals{i});
                end
            else
                for i=1:3
                    obj.Coeffs(i) = vals(i);
                end
            end
        end
%%
        function params = getRunProperties(obj)
            %return properties (open and Transient) as a struct
            mc = metaclass(obj);
			mp = mc.PropertyList;
            ms = mp(1).DefiningClass.SuperclassList;
            scnames = getSuperclassNames(obj,ms);
            count=1;
            propnames = {};
            for k=1:length(mp)
                %remove hidden, transient and constant properties
                idx = mp(k).Hidden;
                %remove superclass properties
                idx = idx + any(strcmp(scnames, mp(k).Name));
                if idx<1
                    propnames{count,1} = mp(k).Name; %#ok<AGROW>
                    count = count+1;
                end
            end
            
            for i=1:length(propnames)
                params.(propnames{i}) = obj.(propnames{i});
            end
            
        end
%% ------------------------------------------------------------------------
% run BMVmodel simulation
%--------------------------------------------------------------------------                
function output =runBMVsimulation(obj,mnvals,mxmn,pvobj)
    %input wave conditions
    varnames = fieldnames(mnvals);
    Hsi = mnvals.(varnames{1}); %mean wave height over averaging interval 
    Tp = mnvals.(varnames{2});  %mean wave period over averaging interval  

   %get "reference" profile for conditions that give lowest closure depth
    if obj.ClosureOption==1
        [maxTp,imx] = max(Tp);
        maxHs = Hsi(imx);
    else
        [maxHs,imx] = max(Hsi);      %maximum wave condition in time series
        maxTp = Tp(imx);
    end
%     K0mx = bmvFittingCoeffs(obj,maxHs,maxTp);
   [obj.Xref,obj.Zref] = referenceProfile(obj,maxHs,maxTp,[]); 
   
   if any(~isreal(obj.Xref))
       msg = sprintf('One or more points of the Minimum Profile are imaginary\nCheck parameter settings');
       warndlg(msg);
       output = [];
       return;
   end
    %
    for i=1:length(Hsi)
        if Hsi(i)>0        
            [X(:,i),Z(:,i)] = sampleProfile(obj,Hsi(i),Tp(i),[]); 
            if mopt==1
                obj.Xref = X(:,i);
                obj.Zref = Z(:,i);
            end
        else
            X(:,i) = NaN; Z(:,i) = NaN;
        end
    end     
    pid = pvobj.mtsc{1}.TimeInfo.UserData;
    titletxt = sprintf('Model profiles for %s',pid);
    legtxt = cellstr(mnvals.date);
    plotEqProfiles(X,Z,Xa,Za,legtxt,omj.Zhw,obj.Zlw,mxmn,titletxt);

    %get profiles within the defined box
    [X,Z,xrange,zrange] = getSampledDataPoints(pvobj,X',Z',mxmn);    
    %get the dimensional volumes and centroids (isND=false)
    [V,X1,Z1] = getVolumeCentroids(pvobj,X,Z,xrange,zrange,false);
    %get the non-dimensional volumens and centroids (isND=true)
    [m0,x1,z1] = getVolumeCentroids(pvobj,X,Z,xrange,zrange,true); 

    output = {V,m0,x1,z1,xrange./zrange,X1,Z1,xrange,zrange}; 
end
%% ------------------------------------------------------------------------
% BMV model functions
%--------------------------------------------------------------------------             
        function [X,Z] = referenceProfile(obj,H,T,K)
            %For the reference profile the distance to high water has to be
            %defined so that the offset between the shoreface and backshore
            %profiles can be established
            if nargin<4 || isempty(K)
                K = bmvFittingCoeffs(obj,H,T);
            end
            %unpack crest position             
            if length(obj.XZc)>1
                xc = obj.XZc(1); obj.Zc = obj.XZc(2);
            else
                xc = 0; obj.Zc = obj.XZc;
            end
            
            [X,Z] =shorefaceProfile(obj,H,T,K);
            X = (X+obj.xHWmn).*(Z<obj.Zhw);
            Xbs = backshoreProfile(obj,H,T,K,Z,obj.xHWmn);
            X = X+Xbs+xc;   %combine the shoreface and backshore profiles
            X = [X;0];      %add point for start of profile
            Z = [Z;obj.Zc+0.1];
        end
%%
        function [X,Z] = sampleProfile(obj,H,T,K)
            %
            if nargin<4 || isempty(K)
                K = bmvFittingCoeffs(obj,H,T);
            end
            %unpack crest position             
            if length(obj.XZc)>1
                xc = obj.XZc(1); obj.Zc = obj.XZc(2);
            else
                xc = 0; obj.Zc = obj.XZc;
            end
            
            [X,Z] = shorefaceProfile(obj,H,T,K);
                
            [ha,hr] = getProfileDepths(obj,H,T,K);
            Za = obj.Zlw - ha;   %closure elevation for each wave condition
            %find adjustment needed to align profiles
            xa = max(X);
            xa_ref = max(obj.Xref);
            Za_ref = min(obj.Zref);
            if xa_ref>xa
                %the new closure distance is less than reference profile 
                xa_icept = interp1(obj.Zref,obj.Xref,Za);
%                 if xa_icept-xc>xa
                    xhw = xa_icept-xa-xc; %offset at Za equates to distance to HW
%                 else
%                     xhw = xa_icept-xa-xc;
%                     
% %                     xhw = xa-xa_icept-xc;
%                 end
            else
                %the new closure distance is greater than reference profile
                xa_icept = interp1(Z,X,Za_ref);
                xhw = xa_ref-xa_icept-xc; %offset at Za_ref equates to distance to HW
            end
 
            X = (X+xhw).*(Z<obj.Zhw); %move x=0 at HW to xhw for shoreface profile 
            Xbs = backshoreProfile(obj,H,T,K,Z,xhw);  %X co-ords for the backshore profle
            X = X+Xbs+xc;  %combine the shoreface and backshore profiles
            X = [X;0];     %add point for start of profile
            Z = [Z;obj.Zc+0.1];  
        end
%%
        function [K,dfv] = bmvFittingCoeffs(obj,H,T)
            %use the settling velocity to calculate the dimensionless fall velocity
            %get system variables
            gr  = obj.g;                   %gravity (m/s^2)
            rhw = obj.rhow;                %water density (kg/m3)
            rhs = obj.rhos;                %sediment density (kg/m3)
            vsc = obj.visc;                %viscosity (m2/s)
            D50 = obj.d50;                 %sediment grain size (m)
            ws  = settling_velocity(D50,gr,rhw,rhs,vsc);%sediment fall velocity (m/s)
            dfv = H./ws./T;                %dimensionless fall velocity (-)
            eqC = eqProfileCoeffs(obj,dfv);%bmv shoreface coefficients
            K = cat(2,obj.Coeffs,eqC);     %backshore + shoreface fitting coefficients
        end
%%
        function [ha,hr] = getProfileDepths(obj,H,T,K)
            %select function to define closure depth, ha
            if obj.ClosureOption==1
                %the shallow water limit is approximated by (h/L=0.09)
                ha = 0.008*obj.g*T^2;
                hr = K(4)*H;
            elseif obj.ClosureOption==2
                %Bernabeu et al use 3Hs
                ha = 3*H;  
                hr = K(4)*H;
            else
                %Hallermeier (1981) 
                mH = obj.meanWaves.mH; 
                sH = obj.meanWaves.sH; 
                mT = obj.meanWaves.mT;
                [ha,hr] = hallermeier_zones(mH,mT,sH,obj.d50,obj.rhow,obj.rhos);
            end       
        end
%%
        function [X,Z] = shorefaceProfile(p,H,T,K)
            %compute the coordinates of the shoreface profile: X,Z
            %find X for Z points (all mOD) at 0.05m intervals
            % p - is an instance of the BMVmodel with defined properties
            % H - Wave height
            % T - wave period
            % K  - fitting coefficients vector

            % hr - breaking wave depth (applied at low water)
            % ha - profile closure depth
            [ha,hr] = getProfileDepths(p,H,T,K); 
            Zr = p.Zlw - hr;   %surf-shoal transition elevation
            Za = p.Zlw - ha;   %closure elevation for each wave condition
            dr = p.Zhw-Zr;     %depth to breaking at high water
            
            %compute X coordinate over range of Z
            pint = 0.05;           %vertical interval for constructed profile (m)
            Z = (Za:pint:p.Zc+pint).';
            
            A = K(5); B = K(6); C = K(7); D = K(8); 
            Xsurf  = (((p.Zhw-Z)/A).^1.5+(B/A^1.5)*(p.Zhw-Z).^3).*(Z<=p.Zhw & Z>=Zr);
            Xo = ((dr)/A)^1.5-(hr/C)^1.5+(B/A^1.5)*(dr)^3-(D/C^1.5)*(hr)^3;  
            Xshoal = (Xo+((p.Zlw-Z)/C).^1.5+(D/C^1.5)*(p.Zlw-Z).^3).*(Z<Zr);
            X = Xsurf+Xshoal;            
        end
%%
        function Xbs = backshoreProfile(p,H,T,K,Z,xhw)
            %compute the coordinates of the backshore profile
            % p - is an instance of the BMVmodel with defined properties
            % H - wave height vector (m)
            % T - wave period (s)
            % K - fitting coefficients vector (use K(3) - dune berm factor)
            % Z - vector of vertical elevation (mOD)
            % xhw - X distance to High water

            Ru = K(2)*runup(p.ubs,H,T);     %wave runup
            Zu = Ru+p.Zhw;                  %upper runup elevation
            Zd = 2*Zu;                      %upper bound of beach profile variation
            H0 = shoaling(H,T,1.28*H,100);  %deep water wave height
            rbs = runup_slope(H0,T,p.ubs,K(1),p.d50,p.S,p.phi,p.g,p.visc);
            if rbs<p.ds/2    %constrain range of run-up slope (Larson,2004 uses 1:1)
                rbs = p.ds/2;
            elseif rbs>100
                rbs = 100;
            end
                
            %variables are defined in terms 2 or 3 characters. The first charater
            %refers to teh x or z dimension, the second to the elevation and tthe
            %third to the feature. Lower case x are distances upper case X are
            %vector arrays for the defined elevations Z. Backshore properties:
            xdd = (p.Zc-Zd)*p.ds;        %distance to dune toe at Zd 
            xul = xdd+(Zd-Zu)*p.ds/2;    %distance to toe of limiting back-berm slope at Zu
            xur = xhw-(Zu-p.Zhw)*rbs;    %X distance to end of runup berm
            dsc = xhw/(p.Zc-p.Zhw);      %slope from crest to high water
            bkc = (xur-xdd)/(Zd-Zu);     %slope from dune toe to top of runup slope
            rbc = (xhw-xul)/(Zu-p.Zhw);  %slope from limiting toe at Zu to high water
            if rbc<rbs 
                %no room for runup slope, steepen rbs to rbc, get updated properties
                xur = xhw-(Zu-p.Zhw)*rbc;%X distance to end of runup berm
                bkc = (xur-xdd)/(Zd-Zu); %slope from dune toe to top of runup slope        
            end

            Rhc = Z>=p.Zhw;              %Z range from high water to crest
            Rhu = Z<Zu & Z>=p.Zhw;       %Z range from high water to top of runup
            Rud = Z<Zd & Z>=Zu;          %Z range from top of runup to dune toe
            Rdc = Z>=Zd;                 %Z range from dune toe to crest

            %identify case and assign slope, range and offset for the 3 parts of
            %the backkshore profile   
            blank = zeros(size(Z));
            drange = Rdc; brange = Rud; rrange = Rhu; %default range settings
            dslope = p.ds; bslope = 0; rslope = 0;    %default slope settings
            if dsc<p.ds/2
                %backshore is at maximum slope (ds/2) from high water to crest
                dslope = p.ds/2;
                drange = Rhc; brange = blank; rrange = blank;
            elseif rbc<p.ds/2
                %backshore is between dune slope and limiting slope (ds/2)
                dslope = dsc;
                drange = Rhc; brange = blank; rrange = blank;
            elseif rbc<rbs
                %dune slope ds back-berm slope ds/2 and runup slope is oversteep
                bslope = p.ds/2;
                rslope = rbc;
                xur = xul;
            elseif bkc<p.ds        
                %dune slope ds, back-berm slope between ds/2 and ds, runup slope rbs
                bslope = bkc;
                rslope = rbs;
            elseif bkc*K(3)<p.ds
                %dune slope ds, back-berm slope ds but limited berm, runup slope rbs
                bslope = p.ds;
                rslope = rbs;
            else
                %dune slope ds, back-berm slope bks with berm, runup slope rbs
                bslope = bkc*K(3);
                rslope = rbs;
            end

            Xdune  = (p.Zc-Z)*dslope.*drange;           %dune profile
            Xback  = (xdd+(Zd-Z)*bslope).*brange;       %back-berm profile
            Xrunup = (xur+(Zu-Z).*rslope).*rrange;      %runup profile      
            Xbs = Xdune+Xback+Xrunup;                   %sum components

            xbs_hw = max(Xbs);
            if xbs_hw>xhw
                %dune has eroded and profile needs to be retreated to high water
                Xbs = (Xbs-xbs_hw+xhw).*(Z>p.Zhw);
            end
        end
%%
        function K = eqProfileCoeffs(obj,dfv)
            %get the empirical coefficients of Bernabeu et al
            % dfv - dimensionless fall velocity
            % option - slection of coefficients to use (2 different papers)
            switch obj.BMVcoeffOption
                case 1
                %coeffs in Bernabeu, ECSS, 2003
                A = 0.13-0.01*dfv;
                B = 0.005+0.26*exp(-0.75*dfv);
                C = 0.11+0.025*dfv;
                D = 0.006+0.1*exp(-0.73*dfv);
                case 2
                %coeffs in Bernabeu, MG, 2003
                A = 0.21-0.02*dfv;
                B = 0.89*exp(-1.24*dfv);
                C = 0.06+0.04*dfv;
                D = 0.22*exp(-0.83*dfv);
            end
            K = [A,B,C,D];
        end
%%
        function K = setFitCoeffs(obj,dfv)
            %reset fit coefficients based on their dependency on dimensionless fall
            %velocity, dfv
            nc = 7;
            props = {'Kc1','Kc2','Kc3','Kc4','Kc5','Kc6','Kc7'};  
            K = zeros(nc,1);
            for i=1:nc
                ab = obj.(props{i});
                if length(ab)==1, ab(2)=0; end %handle constant values
                K(i) = ab(1)+ab(2)*dfv;
                if K(i)<0, K(i) = 0; end
            end
        end
%%
        function tabShoreProfile(obj,src,X,Z,titletxt,XZa,XZr,mnwv)
            %plot profile for selected wave condition
            % obj - instance of BMVmodel
            % src - handle for plot (tab or figure)
            % X,Z - shore profile coordinates
            % titletxt - plot title
            % XZa - shoal closure depth coordinates (optional)
            % XZr - surf limit depth coordinates (optional)
            % mnwv - [mean wave height, std wave height, mean wave period]
            ha = findobj(src.Children,'-depth',1);
            delete(ha)
            if nargin<5
                XZa = []; XZr = []; mnwv = [];
            elseif nargin<6
                XZr = []; mnwv = [];
            elseif nargin<7
                mnwv = [];
            end
            ax = axes(src,'Tag','ShoreProfile');    
            hold on     
            plot(X,Z,'-k');
            if ~isempty(XZa)
                plot(ax,XZa(1),XZa(2),'or');
                txtstr = sprintf('  Za=%.1f (mOD)',XZa(2));
                text(ax,XZa(1),XZa(2),txtstr,'FontSize',6)
            end
            if ~isempty(XZr)
                plot(ax,XZr(1),XZr(2),'og');
                txtstr = sprintf('  Zr=%.1f (mOD)',XZr(2));
                text(ax,XZr(1),XZr(2),txtstr,'FontSize',6)
            end
            if length(obj.XZc)>1
                plot(ax,obj.XZc(1),obj.XZc(2),'om');
                txtstr = sprintf('  Zc=%.1f (mOD)',obj.XZc(2));
                text(ax,obj.XZc(1),obj.XZc(2),txtstr,'FontSize',6)
            end              
            %
            xhw = ceil(max(X));
            xlw = floor(min(X));
            ztop = ceil(max(Z));
            
            plot(ax,[xlw,xhw],[obj.Zhw,obj.Zhw],'--b','DisplayName','High water')
            text(ax,xhw,obj.Zhw,'HWL','FontSize',6)
            plot(ax,[xlw,xhw],[obj.Zlw,obj.Zlw],'--b','DisplayName','Low water')
            text(ax,xhw,obj.Zlw,'LWL','FontSize',6)
            hold off 
            hl = legend('Profile','Shoal limit','Surf limit');
            xlabel('Distance (m)')
            ylabel('Elevation (mOD)')
            title(ax,titletxt)
            %add mean wave values if included
            if ~isempty(mnwv)
                txtstr = sprintf('meahHs=%.1f\nstdHs=%.1f\nmeanTp=%.1f',...
                                            mnwv(1),mnwv(2),mnwv(3));
                ht = text(ax,xhw/2,ztop,txtstr,'FontSize',8);
                ht.Units = 'normalized';
                ht.Position(2) = hl.Position(2)+ht.Extent(4);
            end
        end
%%
        function plotEqProfiles(obj,X,Z,legtxt,titletxt,XZa,XZr,mxmn)
            %plot set of profiles
            % X,Z - array of profiles
            % legtxt - cell array of legend text for each profile
            % titletxt - plot title
            % XZa - shoal closure depth coordinates (optional)
            % XZr - surf limit depth coordinates (optional)
            % mxmn - struct of max and min sampling window (optional)
            if nargin<6
                XZa = []; XZr = []; mxmn = [];
            elseif nargin<7
                XZr = []; mxmn = [];
            elseif nargin<8
                mxmn = [];
            end
            hf = figure('Name','Equilibrium Profiles','Tag','PlotFig','Units','normalized');
            %move figure to top right
            hf.Position(1) = 1-hf.Position(3)-0.01;
            hf.Position(2) = 1-hf.Position(4)-0.12;
            axes;    
            hold on     
            for i=1:length(legtxt)
                hl = plot(X(:,i),Z(:,i));
                if ~isempty(XZa)
                    hp = plot(XZa(i,1),XZa(i,2),'o','Color',hl.Color);
                    set(get(get(hp,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude point from legend
                end
                if ~isempty(XZr)
                    hp = plot(XZr(i,1),XZr(i,2),'.','Color',hl.Color,'MarkerSize',8);
                    set(get(get(hp,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude point from legend
                end
            end
            %
            xhw = ceil(max(max(X,[],2)));
            xlw = floor(min(min(X,[],2)));
            plot([xlw,xhw],[obj.Zhw,obj.Zhw],'--b','DisplayName','High water')
            text(xhw,obj.Zhw,'HWL','FontSize',6)
            plot([xlw,xhw],[obj.Zlw,obj.Zlw],'--b','DisplayName','Low water')
            text(xhw,obj.Zlw,'LWL','FontSize',6)
            if ~isempty(mxmn)
                plot([mxmn.xmin,mxmn.xmin,mxmn.xmax,mxmn.xmax,mxmn.xmin],...
                    [mxmn.zmin,mxmn.zmax,mxmn.zmax,mxmn.zmin,mxmn.zmin],'--r',...
                    'DisplayName','Sampling window');
            end
            hold off 
            if length(legtxt)<20
                legend(legtxt,'Location','southwest');
            end
            title(titletxt)
        end
    end
    
end