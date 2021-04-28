function [K,dfv,iprof,obj] = simBMVfitting(obj,mndst,prdst,mxmn)
%
%-------function help------------------------------------------------------
% NAME
%   simMVfitting.m
% PURPOSE
%   Modified version of the 2-section equilibirum profile model proposed by
%   Bernabeu et al, 2003, extended to include back-beach
%   BMVfitting generates a fit for individual profiles or a times series of 
%   the equilibrium profile fit parameters K1-K3 and optionally A-D.
% USAGE
%   [K,dfv] = Sim_BMVfitting(mnvals,prtsc,mxmn,obj)
%INPUTS
%   obj - BMVmodel properties from Sim_BMV_Input and system variables
%       obj.SedSphericity   %sphericity of sediment
%       obj.BeachPorosity   %beach porosity 
%       obj.GrainSize       %sediment grain size (m)
%       obj.UpperDuneSlope  %upper dune slope (1:ds)
%       obj.UpperBeachSlope %typical/average beach slope (1:bs)
%       obj.BeachCrestLevel %crest elevation (mOD)
%       obj.HWlevel         %high water level (mOD)
%       obj.LWlevel         %low water level (mOD)
%       obj.xHWmn           %minimum distance to high water
%     transient properties of BMVmodel
%       obj.K(1)            %setup coefficient (model fit parameter)
%       obj.K(2)            %runup coefficient (model fit parameter)
%       obj.K(3)            %dune berm factor  (model fit parameter)
%       obj.WaveAvInt       %Wave averaging interval (days) - used prior
%                               to setup input wave data  
%       obj.BMVcoeffOption  %choice of empirical equation set in eqProfileCoeffs
%       obj.AdjustMethod    %Beach adjustment method (1 or 2) - used in BMV model
%       obj.BMVfitOption    %fit using 3 or 7 coefficients
%     system properties added to BMVmodel
%       obj.g               %gravity (m/s^2)   
%       obj.visc            %viscosity (m2/s)
%       obj.rhow            %water density (kg/m3)
%       obj.rhos            %sediment density (kg/m3)
%   mndst - dstable of means Hs and Tp over wave averaging interval, parmas.WaveAvInt
%   prtsc - dstable for selected profle
%   mxmn - 
%OUTPUTS
%   K - profile fit paramaters
%   dfv - dimensionless fall velocity
%   iprof - subset of profiles selected for parameter estimation
% NOTES
%   see Bernabeu et al, ECSS, 2003; Bernabeu, MG, 2003
%   and Reis A H and Gama C, Geomorphology, 2010
% SEE ALSO
%   Sim_BMVinput.m, Sim_BMVmodel.m, CT_Simulation.m

%
% Author: Ian Townend
% CoastalSEA (c)June 2020
%--------------------------------------------------------------------------
%

    %input wave conditions
    % varnames = mndst.VariableNames;
    Hsi = mndst.Hsi; %mean wave height over averaging interval 
    Tp = mndst.Tp;  %mean wave period over averaging interval  

    if obj.ClosureOption==1
        [maxTp,imx] = max(Tp);
        maxHs = Hsi(imx);
    else
        [maxHs,imx] = max(Hsi);      %maximum wave condition in time series
        maxTp = Tp(imx);
    end
%     K0mx = bmvFittingCoeffs(obj,maxHs,maxTp);
    [obj.Xref,obj.Zref] = referenceProfile(obj,maxHs,maxTp,[]);

    obj.ismodel = false;  %fitting function call
%     ncoeff  = obj.BMVfitOption;    %fit using 3, 4, 6 or 8 coefficients
%     if fopt==1, ncoeff=3; else, ncoeff=7; end
    


    %get selected beach profile data   
    t = prdst.RowNames;
    x = prdst.Chainage; 
    z = prdst.Elevation;
    pid = prdst.Description;
    
    imxdate = datestr(t(imx),'dd-mmm-yyyy');
    cfProfiles(x(imx,:),z(imx,:),obj.Xref,obj.Zref,obj.Zhw,obj.Zlw,mxmn,imxdate,'test');
    %get user to select a profile (or all profiles) to fit
    ncoeff = 8;  %number of coefficients in K
    ok = 1;   K = [];   dfv = []; iprof = [];
    tlist = datestr(t);
    fprintf('Limiting Hs = %.2f and Tp = %.2f on %s\n',maxHs,maxTp,imxdate)
    %
    while ok>0
        [h_dlg,ok] = listdlg('Name','Plot profile', ...
            'PromptString','Select survey date(s)', ...
            'SelectionMode','multiple', ...
            'ListString',tlist);
        if ok==0, return; end  %user cancelled
        nrec = length(h_dlg)-sum(isnan(Hsi(h_dlg)));
        if nrec>1
            %analyse profiles and return timeseries of K and dfv
            ok = 0;
            dfv = NaN(nrec,1); iprof = dfv;
            K = NaN(nrec,ncoeff);
            aflag = 1; %flag to indicate multiple profiles being fitted
            ic = 1;    %counter so that any profiles can be selected as set
            hw = waitbar(0, 'Processing.');
            for ip=h_dlg
                date = datestr(t{ip},'dd-mmm-yyyy');
                if ~isnan(Hsi(ip))
                    txt = sprintf('comparison for profile %s surveyed on %s',pid,date);
                    [K(ic,:),dfv(ic,1)] = fitProfile(obj,Hsi(ip),Tp(ip),...
                                mxmn,x(ip,:)',z(ip,:)',txt,aflag);
                    iprof(ic,1) = ip;
                    ic = ic+1;
                end
                fprintf('%g: Profile %s surveyed on %s\n',ip,pid,date);
                waitbar(ip/length(h_dlg))
            end
            close(hw)
        else
            %select profiles individually
            if ok==0,  return; end
            K = 999; 
            aflag = 0;   %flag to indicate single profile being fitted
            date = datestr(t(h_dlg),'dd-mmm-yyyy');
            while K==999
                xi = (x(h_dlg,:)');
                zi = (z(h_dlg,:)');
                txt = sprintf('Model comparison for profile %s surveyed on %s',pid,date);
                [K,dfv] = fitProfile(obj,Hsi(h_dlg),Tp(h_dlg),mxmn,...
                                                xi,zi,txt,aflag);
                iprof = h_dlg;
                if isempty(K) 
                    return; 
                elseif K==999
                    idx = 1:9;
                    obj = editPropertySubset(obj,idx);
                    obj = setCoeffProps(obj);
                    obj = setRunProps(obj);
                else
                    %display updated estimates of coefficients
                    fopt  = obj.BMVfitOption;    %fit using 3, 4, 6 or 8 coefficients
                    if fopt==3
                        msgtxt = sprintf('Updated estimates for profile %s\nSurveyed on %s:\nK1: %g\nK2: %g\nK3: %g',...
                           pid,date,K(1),K(2),K(3));
                    elseif fopt==4
                        msgtxt = sprintf('Updated estimates for profile %s\nSurveyed on %s:\nK1: %g\nK2: %g\nK3: %g\nK4: %g\n',...
                           pid,date,K(1),K(2),K(3),K(4));
                    elseif fopt==6  
                        msgtxt = sprintf('Updated estimates for profile %s\nSurveyed on %s:\nK1: %g\nK2: %g\nK3: %g\nK4: %g\nK5: %g\nK6: %g\n',...
                           pid,date,K(1),K(2),K(3),K(4),K(5),K(6));
                    else
                        msgtxt = sprintf('Updated estimates for profile %s\nSurveyed on %s:\nK1: %g\nK2: %g\nK3: %g\nK4: %g\nK5: %g\nK6: %g\nK7: %g\nK8: %g\n',...
                           pid,date,K(1),K(2),K(3),K(4),K(5),K(6),K(7),K(8));
                    end
                    h_m = msgbox(msgtxt);
                    h_b = findobj(h_m,'String','OK');
                    but_pos = h_b.Position;
                    but_pos(1) = but_pos(1)+ 1.2*but_pos(3);
                    %Create push button to copy data to clipboard
                    uicontrol('Parent',h_m,...
                        'Style','pushbutton',...
                        'String', sprintf('Copy'),...
                        'Units',h_b.Units, ...
                        'Position', but_pos, ...
                        'UserData',msgtxt, ...
                        'Callback',@(src,evdata)mat2clip(src.UserData));
                    uiwait(h_m);
                end
            end
        end
    end    
end
%%
function [K,dfv] = fitProfile(obj,Hs,Tp,mxmn,x,z,txt,aflag)
    %fit the selected measured profile to the BMV profile

    %remove NaNs from profile
    prof = table(x,z);
    prof = rmmissing(prof,1); %remove rows with NaNs
    x = table2array(prof(:,1));
    z = table2array(prof(:,2));
    [x,z] = rangeProfile(x,z,mxmn,2);      
   
    %default values of K1-K3, gamma and A-D: calls function eqProfileCoeffs in bmvprofile 
    [K0,dfv] = bmvFittingCoeffs(obj,Hs,Tp);
    fc = 5;
    %constrained optimisation upper and lower bounds
    fopt= obj.BMVfitOption;      %fit using 3, 4, 6 or 8 coefficients
    if fopt==3   %fit using just the shape coefficients
        ub = [0.5, 1.5,  1, 1,1,1,1,1];
        lb = [0.1, 0.5,  0, 1,1,1,1,1];   
    elseif fopt==4
        ub = [0.5, 1.5,  1, 1.5,  1,1,1,1];
        lb = [0.1, 0.5,  0, 0.5,  1,1,1,1];   
    elseif fopt==6%fit using shape coefficients and surf zone dissipation coefficients
        ub = [0.5, 1.5,  1, 1.5,  fc*K0(5),  fc*K0(6), 1,1];
        lb = [0.1, 0.5,  0, 0.5,  eps,       0,        1,1];  
    else         %fit using shape coefficients and all dissipation coefficients
        ub = [0.5, 1.5,  1, 1.5,  fc*K0(5),  fc*K0(6), fc*K0(7), fc*K0(8)];
        lb = [0.1, 0.5,  0, 0.5,  eps,       0,        eps,      0];       
    end
    
    options = optimoptions('lsqcurvefit','Display','off');
    %call for testing and debugging----------------------------------------
    %options = optimoptions('lsqcurvefit','Display','iter-detailed',...
    %                       'Diagnostics','on','FunctionTolerance',1e-6,...
    %                       'StepTolerance',1e-6,'OptimalityTolerance',1e-6);
    %----------------------------------------------------------------------
    %get minimum "reference" profile
    [maxH,imx] = max(Hs);  %maximum wave condition in time series
    [obj.Xref,obj.Zref] = referenceProfile(obj,maxH,Tp(imx),K0);
    %now use this to fit against
    interpProfile = @(K,z) fitSampleProfile(obj,Hs,Tp,K,z);
    
    if aflag
        K = lsqcurvefit(interpProfile,K0,z,x,lb,ub,options);  
    else
        X = interpProfile(K0,z);    
        legtxt = {'Observed','BMV model'};        
        hf = cfProfiles(x,z,X,z,obj.Zhw,obj.Zlw,mxmn,legtxt,txt);

        ansQ = questdlg('Accept initial fit and run regression?',...
                        'BMV fitting','Yes','Adjust','Quit','Yes');
        if strcmp(ansQ,'Yes')
            K = lsqcurvefit(interpProfile,K0,z,x,lb,ub,options);
            X = interpProfile(K,z);           
            cfProfiles(x,z,X,z,obj.Zhw,obj.Zlw,mxmn,legtxt,txt,hf);
        elseif strcmp(ansQ,'Adjust')
            delete(hf)
            K = 999;  return;
        else
            delete(hf)
            K = [];   return;
        end
    end
end
%%
        function Xz = fitSampleProfile(obj,H,T,K,z)
            %called from Sim_BMVfitting to get a sampleProfile and
            %interepolate onto measured profile elevation intervals, z
            K0 = bmvFittingCoeffs(obj,H,T);
            fopt = obj.BMVfitOption;   %fit using 3, 4, 6 or 8 coefficients
            if fopt==3
                K(4:8) = K0(4:8);  %use only backshore coefficients
            elseif fopt==4
                K(5:8) = K0(5:8);  %use backshore and breaker index
            elseif fopt==6
                K(7:8) = K0(7:8);  %use backshore, breaker + surf coefficents
            elseif fopt==8
                                   %use all K values (including shoaling)
            else
                warndlg('Invalid number of fit parameters')
                Xz = []; return;
            end

            [X,Z] = sampleProfile(obj,H,T,K);
            Xz = interp1(Z,X,z);
            Xz(isnan(Xz)) = 0;
        end
%%
function [xi,zi] = rangeProfile(vx,vz,mxmn,flag)
    %find points at each end of profile nearest to xmin and zmin
    %flag 0=values from x=0,z=0; 1= normalised values; 2=values to OD
    xmin = mxmn.xmin;                %profile sampling window x minimum
    zmin = mxmn.zmin;                %profile sampling window z minimum
    
    [~,Ix] = min(abs(vx-xmin));
    [~,Iz] = min(abs(vz-zmin));
    %select data points above minimum values
    idv = find(vx>=vx(Ix) & vz>=vz(Iz));
    xi = vx(idv)-vx(Ix);
    zi = vz(idv)-vz(Iz);
    %zi = zi-zmin;           %make elevation heights above zmin
    %handle nans in data
    idx = isnan(xi);
    xi(idx) = 0;
    zi(idx) = 0;
    
    switch flag
        case 1   %option to normalise data to range of x and z
            maxXrange = max(vx(Iz)-vx(Ix));
            maxZrange = max(vz(Ix)-vz(Iz));
            xi = xi/maxXrange;
            zi = zi/maxZrange;
        case 2   %option to return values relative to Datum (OD)
            xi = xi+xmin;
            zi = zi+zmin;
    end
end
%%
function hf = cfProfiles(X1,Z1,X2,Z2,Zhw,Zlw,mxmn,legtxt,titletxt,hf)
    %plot profiles
    if nargin<10
        hf = figure('Name','Equilibrium Profiles','Tag','PlotFig','Units','normalized');
        %move figure to top right
        hf.Position(1) = 1-hf.Position(3)-0.01;
        hf.Position(2) = 1-hf.Position(4)-0.12;
        axes; 
        hold on
        plot(X1,Z1,'DisplayName','Surveyed profile');
        plot(X2,Z2,'DisplayName','Initial profile');
        %add high and low water lines
        xhw = ceil(max(max(X1),max(X2)));
        xlw = floor(min(min(X1),min(X2)));
        plot([xlw,xhw],[Zhw,Zhw],'--b','DisplayName','High water')
        plot([xlw,xhw],[Zlw,Zlw],'--b','DisplayName','Low water')
        %plot survey sampling window
        plot([mxmn.xmin,mxmn.xmin,mxmn.xmax,mxmn.xmax,mxmn.xmin],...
         [mxmn.zmin,mxmn.zmax,mxmn.zmax,mxmn.zmin,mxmn.zmin],'--r',...
         'DisplayName','Sampling window');
        if length(legtxt)<20
            legend(legtxt,'Location','northeast');
        end
        title(titletxt);
    else
        hold on
        plot(X2,Z2,'DisplayName','Fitted profile');
    end 
    hold off 
end