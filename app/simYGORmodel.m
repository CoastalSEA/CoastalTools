function [tt,ss,mtxt] = simYGORmodel(wv,sh,wvb,coeffs,isfcast)     
%
%-------function help------------------------------------------------------
% NAME
%   Sim_YGORmodel.m
% PURPOSE
%   Function to clean data and call optimization algorithm to find
%   the best fit for the five free parameters in the Yates et al model

% USAGE
%   [tt,ss,mtxt] = Sim_YGORmodel(wv,b1,wvb,labels,coeffs,isfcast)   
% INPUTS
%   wv  - wave parameter (eg energy, etc) as a dstable
%   sh  - beach profile parameter (eg shoreline position) as a dstable
%   wvb - averaged wave parameter for preceding interval of beach profile
%   coeffs - fit coefficients to be used in forecast
%   isfcast - flag for forecast rather than model run
% OUTPUT
%   When called to fit data
%     tt - the fitted coefficients
%     ss - estimates of the root means square error (rmse)
%     mtxt - not used (returns [])
%   when called to make a simulation
%     tt - time
%     ss - estimate of variable using defined model fit coefficients
%     mtxt - metatext that summarises the simulation
% NOTES
%   see Yates et al, JGR, 2009 and Villamarin, MSc thesis, UoS, 2017
%   Requires the Matlab Global Optimization Toolbox
% SEE ALSO
%   Sim_YGORinput.m, CT_Simulation.m
%
% Author: Ian Townend
% CoastalSEA (c)June 2017
%--------------------------------------------------------------------------
%

    %remove NaNs for shoreline data and energy between profiles
    [ss,ts,Es,Et,tt,idx] = cleandata(sh,wvb,wv,isfcast);
    labels = [sh.VariableDescriptions(1),wv.VariableDescriptions(1)];
    %remove any underlying trend in the data
    dts = days(ts-ts(1));
    pfitss = polyfit(dts,ss,1);
    trend = polyval(pfitss,dts);
    ssdt = ss-trend;
    %ssdt = ss;
    %check code------------------------------------------------------------
    %detrend the shoreline position data (this creates a ts about trend line)
    %pfitssdt = polyfit(dts,ssdt,1);
    %figure(1)
    %plot(ts, ss,'-b',ts,polyval(pfitss,dts),'--r')
    %hold on
    %plot(ts,ssdt, '-k',ts,polyval(pfitssdt,dts),'--r')
    %hold off
    %----------------------------------------------------------------------
    
    if isfcast
        [tt,ss,mtxt] = getYGORfcast(ss,ts,Et,tt,coeffs,pfitss,labels);
    else
        mtxt = [];
        [tt,ss] = getYGORfit(ss,ts,Es,Et,tt,ssdt,idx,trend,labels);
    end
end
%%
function [lastC,rmse] = getYGORfit(ss,ts,Es,Et,tt,ssdt,idx,trend,labels)
    %fit shoreline change to the model forcing parameter
    %use values between profiles to get initial guess for a and b
    lastC = []; rmse = []; 
    [p1,p2,dsdt] = dSdt_plot(ssdt,Es,ts,labels);     %#ok<ASGLU>
    if isempty(p2), return; end
    %to use the fit to all the data to inialise a and b use p1 and to use
    %the subset for dsdt->0 defined in dSdt_plot use p2
    %a0 = p1(1); b0 = p1(2);
    a0 = p2(1); b0 = p2(2);
%     if a0<abs(1e-4)   %very small values of 'a' narrow the search too much
%         a0 = sign(a0)*1e-4;
%     end
    [a,b] = adjustAB(a0,b0); %prompt user to accept estimates or adjust them

    %prompt user to select an optimisation algorithm
    AlgList = {'fminsearch','simulanneal','particleswarm','patternsearch'};
    [select, ok] = listdlg('PromptString','Select an algorithm:',...
                'SelectionMode','single','ListSize',[150,65],'ListString',AlgList);
    if ok<1, return; end  %user aborts

    %call function to compute model parameters
    [C, fval] = YGORparameters(ssdt,Et,tt,idx,a,b,dsdt,AlgList{select});
    fprintf('Final estimates: a=%.3g, b=%.3g, C(+)=%.3g, C(-)=%.3g, d=%.3g fval=%.3g\n',...
        C(1),C(2),C(3),C(4),C(5),fval);
    
    %now plot actual and simulated shoreline change
    while ~isempty(C)    
        [tx,~,idt] = intersect(ts,tt,'rows');
        Tt = days(tt-tt(1));           %time as a numeric array from t=0        
        [St, ~] = YGORts(ssdt(1,1),Et,Tt,C);
        smod = St(idt); 
        smod = smod+trend;
        %delete records with NaNs in both datasets first
        %I = ~isnan(ss) & ~isnan(modata); 
        %ss = ss(I); smod = smod(I);
        R = corrcoef(ss,smod);
        R2 = R(1,2)^2;
        rmse =sqrt(sum((ss(:)-smod(:)).^2)/numel(ss));
        
        YGORplot(ts,ss,tx,smod,labels,C,R,R2,rmse);
        lastC = C;
        C = getC(C);        
    end
    %plot goodness of fit
    YGORfitplot(ss,smod,lastC,R2);
    [skill,ok] = setskillparameters();
    if ok>0 && skill.Inc             
        metatxt{1} = 'Data';
        metatxt{2} = 'Model';
        radial_limit = 2;
        taylor_plot(ss,smod,metatxt,'New',radial_limit,skill);
    end
end
%%
function [tt,smod,mtxt] = getYGORfcast(ss,ts,Et,tt,coeffs,pfitss,labels)
    %generate forecast based on defined coefficients
    tf = days(tt-tt(1));           %time as a numeric array from t=0     
    [smod, ~] = YGORts(0,Et,tf,coeffs.C);
    smod = smod+polyval(pfitss,tf);  
    %adjust position to a defined point at offsetdate       
    offsetdate = datetime([coeffs.Dat0,' 00:00:00'],'InputFormat','dd-MM-yyyy HH:mm:ss');
    idt = tt==offsetdate;
    if any(idt)
       offset = smod(idt)-coeffs.Sh0;
    else
       offset = double(idt);
    end
    smod = smod-offset;
    %add startdate to label
    startdate = datestr(tt(1),1);
    labels{2} = sprintf('%s from %s', labels{2}, startdate);
    %plotting adaptation-----------------------------------------------
    %convert to years so that can use "patch" as backdrop (which needs
    %doubles not datetime (used in ShoreCast)
    %startyear = year(tt(1));
    %ts = startyear+years(ts-datetime(startyear,1,1)); 
    %tt = startyear+years(tt-datetime(startyear,1,1)); 
    %------------------------------------------------------------------
    if ts(end)<tt(1)  %if shoreline data does not overlap forecast period
        ts = []; ss = [];
    end
    YGORplot(ts,ss,tt,smod,labels,coeffs.C);
    ctxt = sprintf('a=%.3g, b=%.3g, C(+)=%.3g, C(-)=%.3g, d=%.3g',...
              coeffs.C(1),coeffs.C(2),coeffs.C(3),coeffs.C(4),coeffs.C(5));
    mtxt = sprintf('%s, %s, %s',labels{1},labels{2},ctxt);
end
%%
function YGORplot(ts,ss,tx,smod,labels,C,R,R2,rmse)
    %plot model predictions against the measured data
    if nargin<7
        R2 = -99;
    end
    hf = figure('Name','YGOR model','Tag','PlotFig','Units','normalized');
    %move figure to top right
    hf.Position(1) = 1-hf.Position(3)-0.01;
    hf.Position(2) = 1-hf.Position(4)-0.12;
    plot(ts,ss,'DisplayName','Observed')
    hold on
    plot(tx,smod,':.','DisplayName','Eq(2)') 
    %plot(tx,Ss(idt)+trend,':+','DisplayName','Eq(6)')  
    %legend({'Observed','Model Eq(2)','Model Eq(6)'},'Location','best')
    %Sntrend = polyval(pfitss,Tt);
    %plot(tt,St+Sntrend,':.','DisplayName','Eq(2)')  %full time series        
    hold off
    legtxt1 = sprintf('Observed using %s',labels{1});
    legtxt2 = sprintf('Model using %s',labels{2});
    if isempty(ss)
        legend({legtxt2},'Location','best')
    else
        legend({legtxt1,legtxt2},'Location','best')
    end
    if R2<0
        titletxt = sprintf('Model using: a=%.3g, b=%.3g, C(+)=%.3g, C(-)=%.3g, d=%.3g',...
                            C(1),C(2),C(3),C(4),C(5));
    else
        titletxt = sprintf('Model using: a=%.3g, b=%.3g, C(+)=%.3g, C(-)=%.3g, d=%.3g\nR=%g; R^2=%0.3g; RMSE=%0.3g',...
                            C(1),C(2),C(3),C(4),C(5),R(1,2),R2,rmse);
    end
    title(titletxt);
    xlabel('Time');
    ylabel(labels{1});
end
%%
function YGORfitplot(ss,smod,C,R2)
    %plot measured v modelled with 1:1 and best fit lines
    smax = max(ss); smin = min(ss);
    smmx = max(smod); smmn = min(smod);
    minpt = min(smin,smmn); maxpt = max(smax,smmx);
    
    hf = figure('Name','YGOR model','Tag','PlotFig','Units','normalized');
    %move figure to top right
    hf.Position(1) = 1-hf.Position(3)-0.01;
    hf.Position(2) = 1-hf.Position(4)-0.12;
    dotsize = 5;
    scatter(ss,smod,dotsize,'filled')
    hold on
    k1 = lsline;
    plot([minpt,maxpt],[minpt,maxpt],'--k');
    hold off
    %set line properties and get fit coefficients
    set(k1(1),'color','r','LineStyle','-.'); % line properties to be visible
    p1 = polyfit(get(k1(1),'xdata'),get(k1(1),'ydata'),1);  %all data
   
    %add the line coefficients to the plot
    tx1 = 'Fit coefficients:';
    txtstr = sprintf('%s\nSlope: %.3g\nIntercept: %.3g\nR^2: %.3g',...
                         tx1,p1(1),p1(2),R2);
    % position to display in the scatter plot
    xpos=min(ss)-3; ypos=max(smod)-3;
    text(xpos,ypos,txtstr,'FontSize',9) 
    titletxt = sprintf('Model using: a=%.3g, b=%.3g, C(+)=%.3g, C(-)=%.3g, d=%.3g',...
                            C(1),C(2),C(3),C(4),C(5));
    title(titletxt);
    xlabel('Observed shoreline position (m)');
    ylabel('Modelled shohreline position (m)');
    legend('Data','Best fit','1:1','Location','best');
end
%%
function [S,Ss] = YGORts(s0,E,t,C)
%use the fitted parameters to reconstruct timeseries of shoreline change
    Seq = (E-C(2))/C(1);
    S(1,1) = s0+C(5);   Ss(1,1) = s0+C(5); 
    for i=2:length(t)
        Eeq = C(1)*S(i-1,1)+C(2);       %Eeq(S)=aS+b  eq(4)
        delE = E(i)-Eeq;                %delE(S) = E-Eeq(S)  eq(3)
        if delE<0
            Cpm = C(3);                 %accretion rate
        else
            Cpm = C(4);                 %erosion rate
        end
        %shoreline change using solution as per eq(6)-ONLY valid for CONSTANT E
        Ss(i,1) = (Ss(1,1)-Seq(i)).*exp(-C(1)*Cpm*sqrt(E(i))*t(i))+Seq(i);
        %shoreline change using rate of shoreline change eq(2):
        S(i,1) = S(i-1,1)+Cpm*sqrt(E(i))*delE;        
    end
end

%%
function S = YGORsubsample(s0,E,t,C,idx)
    %subsample time series based on index values idx
    St = YGORts(s0,E,t,C);
    S  = St(idx);
end

%%
function [C, fval] = YGORparameters(ss,Et,tt,idx,a,b,dsdt,OptAlg)
    %use derivative free, constrained non-linear optimisation to find the best
    %fit for the five parameters in the YGOR model: a, b, C+, C- and origin
    %offset
    %inputs are beach variable (eg position), energy and optimisation algorithm    
    %ss - cleaned shore variable at times ts
    %Et - 'wave energy' at time tt (full time series)
    %idx - index for profiles in wave timeseries
    %a,b - initial guess fit parameters
    %dsdt - positive and negative rate of shoreline change
    
    Tt = days(tt-tt(1));        %time for waves as a numeric array from t=0
    s0 = ss(1,1);  %initial shoreline position

    funSx = @YGORsubsample;
    funS = @(C) funSx(s0,Et,Tt,C,idx);
    funRMSE = @(C) sqrt(mean(funS(C)-ss).^2); %objective function to solve
    
    Co(1) = a;              %fit parameter for equilibirum slope, a
    Co(2) = b;              %fit parameter for equilibirum intercept, b
    Co(3) = -dsdt(2);       %accretion rate coefficient for deltaE<0, C+
    Co(4) = dsdt(1);        %erosion rate coefficient for deltaE>0, C-
    Co(5) = 0;              %offset for initial position of shoreline
   
    %constrained optimisation upper and lower bounds    
    ub = [-0.5*abs(a),  1.2*b,  -0.1,     -0.1,        1.2*std(ss,'omitnan')];
    lb = [-1.5*abs(a),  0.8*b,   2*Co(3),  2*Co(4),   -1.2*std(ss,'omitnan')]; 
    
    switch OptAlg
        case 'fminsearch'
            %unconstrained optimisation
            options = optimset('MaxIter',2000,...
                'TolFun',1e-5,'TolX',1e-9'); %,...
                %'PlotFcns',{@optimplotx,@optimplotfval});    
            [C, fval] = fminsearch(funRMSE,Co,options); %sensitive to initial values

        case 'simulanneal'    
            %constrained optimisation using simulated annealing
            %options plots the process but also slows the function down
            %- temperature plot is @saplottemperature
            options = optimoptions(@simulannealbnd,...
                    'MaxIterations',2000,...
                    'FunctionTolerance',1e-5,...
                    'ObjectiveLimit',1e-9,...
                    'InitialTemperature',[150,50,200,200,100],...
                    'TemperatureFcn',@temperatureexp,...
                    'PlotFcns',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
            [C, fval] = simulannealbnd(funRMSE,Co,lb,ub,options);
            
        case 'particleswarm'
            %constrained optimisation using particle swarm
            options = optimoptions(@particleswarm,...
                    'PlotFcn',@pswplotbestf);
            [C, fval] = particleswarm(funRMSE,length(lb),lb,ub,options);

        case 'patternsearch'
            %constrained optimisation using pattern search with linear constraint
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            nonlcon = [];
            options = optimoptions('patternsearch','Display','iter','PlotFcn',...
                {@psplotbestf,@psplotmeshsize,@psplotfuncount,@psplotbestx});
            [C, fval] = patternsearch(funRMSE,Co,A,b,Aeq,beq,lb,ub,nonlcon,options);
    end
end
%%
function [slope,intercept,Rsq] = YGORslope(x,y)
    % see Matlab documentation for Linear Regressioin > Simple Linear Regression
    % Simple linear regression considers only one independent variable using the
    % relation
    %  y = a.x + b, where b is the y-intercept and a is the slope    
    X = [ones(length(x),1) x];
    a_b = X\y;
    yCalc = X*a_b;
    Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);
    slope = a_b(2);
    intercept = a_b(1);
end

%%
function C = getC(C)
    %allow user to modify esitmated values of C
    p1 = sprintf('Press "OK" to change values\nPress "Cancel" to use current values\na: ');
    prompt = {p1,'b:','C(+):','C(-):','d:'};
    title = 'Adjust free parameters';
    numlines = 1;
    defaultvalues = {num2str(C(1)),num2str(C(2)),num2str(C(3)),...
        num2str(C(4)),num2str(C(5))};
    C = inputdlg(prompt,title,numlines,defaultvalues);
    if isempty(C), return; end
    C = str2double(C);
end
%%
function ESplot(s,t,Et,tt,ptype) %#ok<DEFNU>
    %plot profile position and full energy (both normalised)
    %ptype defines whether second variable is line or bar plot
    sn = normalise(s);
    Etn = normalise(Et);
    hf = figure('Name','YGOR data','Tag','PlotFig','Units','normalized');
    %move figure to top right
    hf.Position(1) = 1-hf.Position(3)-0.01;
    hf.Position(2) = 1-hf.Position(4)-0.12;
    bar(t,sn,'r','EdgeColor','r','DisplayName','Position');
    hold on
    if ptype
        bar(tt,Etn,2,'EdgeColor','b','FaceColor','none',...
            'LineWidth',0.2,'DisplayName','Energy');
    else
        plot(tt,Etn,'b','DisplayName','Energy');
    end
    xlabel('Time')
    ylavel('Normalised values')
    legend({'Position','Energy'})
    hold off
end

%%
function newvar = normalise(var)
   %takes the mean and std for each column and normalises each column
    mvar = mean(var,'omitnan');
    svar = std(var,'omitnan');
    newvar = (var-mvar)./svar;
end
%%
function [ss,ts,Es,Et,tt,idx] = cleandata(shdst,wvb,wvdst,isfcast)
    %remove NaNs and clip to duration of profile/shoreline data
    tsh = shdst.RowNames;
    sh = shdst.DataTable{:,1};
    twv = wvdst.RowNames;
    wv = wvdst.DataTable{:,1};
    allprofs = table(tsh,wvb,sh);
    allprofs = rmmissing(allprofs,1); %remove rows with NaNs   
    t = table2array(allprofs(:,1));  %time of profile survey
    E = table2array(allprofs(:,2));  %average "energy" over interval prior to profile
    s = table2array(allprofs(:,3));  %observed shoreline position
    I = ~isnan(wv); wv = wv(I); twv = twv(I);
    %crop full wave record to length of profile time series
    t.Format = 'dd/MM/yyyy HH:mm:ss';    
    if isfcast
        Et = wv; tt = twv;
        ss = s; ts = t; Es = E; idx = [];
        idx = isnan(wv);     
    else
        Et = wv(isbetween(twv,t(1),t(end)));  %clipped energy timeseries
        tt = twv(isbetween(twv,t(1),t(end))); %clipped time of eneryg ts
        tt.Format = 'dd/MM/yyyy HH:mm:ss';
        
        %find profile times in wave time series twv=tt(idx)
        [tx,~,idx] = intersect(t,tt,'rows');         
        if sum(isnan(Et))==0
            method = 4;
            switch method
                case 1  %interpolate waves to time of profile                    
                    Es = interp1(tt,Et,t);
                case 2  %use mean or max of waves between profiles
                    for i=1:length(t)-1
                        %Es(i,1) = mean(Et(isbetween(tt,t(i),t(i+1))));
                        Es(i,1) = max(Et((tt>=t(i)& tt<=t(i+1))));
                    end
                    Es = [Es(1);Es]; %assume waves before fist profile = first interval
                case 3 %use wave at prior profiletimestep
                    Es = interp1(tt,Et,t);
                    Es = [Es(1);Es(1:end-1)]; 
                case 4
                    Es = Et(idx);
            end
            ss = s;
            ts = t;
        else            
            %remove times with profiles but no waves ts = tx(ids) ie not index of tt
            [ts,~,ids] = intersect(t,tx,'rows');  
            ss = s(ids);
            Es = E(ids); 
        end
        %ESplot(ss,ts,Et(idx),tx,true);  %plot of shoreline and subset of waves  
        %ESplot(s,t,Et,tt,false);   %plot source data (plots normalised values)  
    end   
end
%%
function [p1,p2,mdxdt] = dSdt_plot(x,y,t,labels)
    %generate scatter plot of data and the fit lines for dx/dt->0
    
    idy = y<-eps | y>eps;                %find index of values ~=0
    x = x(idy); y = y(idy); t = t(idy);  %remove values very close to 0
    
    if length(x)<20                      %not enough data points
        p1 = []; p2 = [];mdxdt = [];
        warndlg('Not enough shoreline data points to produce a meaningful fit')
        return;
    end
    
    %rate of change of x (eg shoreline position)
    dxdt = [0;diff(x)./days(diff(t))]; %forward stepping gradient
    dxmax = max(max(dxdt),abs(min(dxdt)));
    dxmean = mean(abs(dxdt));
    mdxdt(1) = mean(dxdt(dxdt<0));
    mdxdt(2) = mean(dxdt(dxdt>0)); 
    %fit to all the data to establish a lower threshold for y
    [a,b,rsq] = YGORslope(x,y);  %use avearge values for initial guess
    fprintf('Initial estimates: a = %.3g, b = %.3g, rsq = %.3g\n',a,b,rsq)

    %find values that are close to the minimum of dx/dt (eg rate of change of shoreline)
    minfact = 0.05; miny = 0.01;
    idx = find((dxdt<=minfact*dxmean) & (dxdt>=-minfact*dxmean) & y>miny); 
    %if there are not enough points to fit a line find 10 smallest values
    if length(idx)<11
        [~,isx] = sort(abs(dxdt));
        idx = isx(1:10);
    end
    idx = idx(2:end); %remove leading zero used to pad arrary dxdt
    
    %generate plot
    hf = figure('Name','YGOR dSdt','Tag','PlotFig','Units','normalized');
    %move figure to top right
    hf.Position(1) = 1-hf.Position(3)-0.01;
    hf.Position(2) = 1-hf.Position(4)-0.12;
    dotsize = 10;
    scatter(x,y,dotsize,dxdt,'filled');    %all values
    hold on
    scatter(x(idx),y(idx),2*dotsize); %values close to minimum
    hold off
    xlabel(labels{1});
    ylabel(sprintf('Av.%s',labels{2}));
    h = colorbar;
    ind = regexpi(labels{1},'(');
    ltxt = labels{1}(1:ind-2);
    unitxt = labels{1}(ind:end-1);
    h.Label.String = sprintf('%s change rate %s/day)',ltxt,unitxt);
    load 'cmapanomalie.mat'  cmapanomalie %colormap from red to blue with white as 0.
    colormap(cmapanomalie);
    dxscale = dxmax/2;
    caxis([-dxscale dxscale]);
    set(gcf,'color','w');

    %now get only points with dxdt close to zero
    k1 = lsline;
    hold off
    %set line properties and get fit coefficients
    set(k1(1),'color','k','LineStyle','--'); % line properties to be visible
    p1 = polyfit(get(k1(1),'xdata'),get(k1(1),'ydata'),1);  %all data
    set(k1(2),'color','k','LineWidth',1); % line properties to be visible    
    p2 = polyfit(get(k1(2),'xdata'),get(k1(2),'ydata'),1);  %subset dsdt->0
   
    %add the line coefficients to the plot
    mE = mean(y,'omitnan');
    tx1 = 'Fit coefficients: All (Subset)';
    txtstr = sprintf('%s\nSlope: %.3g (%.3g)\nIntercept: %.3g (%.3g)\nMean Energy: %.3g',...
                         tx1,p1(1),p2(1),p1(2),p2(2),mE);
    % position to display in the scatter plot
    xpos=0; ypos=max(y);
    text(xpos,ypos,txtstr,'FontSize',9)
    legend('Data','Data ds/dt->0','Fit to all','Fit to ds/dt->0','Location','northwest');
end
%%
function [a,b] = adjustAB(a0,b0)
    %prompt user to accept estimates of a and or adjust them
    while ~isempty(a0)
        figax = gca;
        figline = findobj(figax,'Type','line');
        figline = findobj(figline,'Tag','userline');
        if ~isempty(figline)
            delete(figline)
        end
        hold on
        xx = figax.XTick;
        yy = b0+a0*xx;
        xx = xx(yy>=0);
        yy = yy(yy>=0);        
        plot(xx,yy,'--r','LineWidth',1,'DisplayName','Selected Fit');
        hold off
        a = a0; b = b0;
        [a0,b0] = getAB(a,b);
    end
end
%%
function [a,b] = getAB(a,b)
    %allow user to modify esitmated values of C
    p1 = sprintf('Press "OK" to change values\nPress "Cancel" to use current values\na: ');
    prompt = {p1,'b:'};
    title = 'Adjust fit parameters';
    dims = [1 40];
    defaultvalues = {num2str(a),num2str(b)};
    C = inputdlg(prompt,title,dims,defaultvalues);
    if isempty(C), a = []; return; end
    a = str2double(C{1});
    b = str2double(C{2});
end