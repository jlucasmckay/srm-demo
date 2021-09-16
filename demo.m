function demo()
% demo()
% 
% demo.m recreates Figure 2A in:
% 
% McKay et al., "Abnormal center of mass control during balance: a new
% biomarker of falls in people with Parkinson’s disease"
% 
% https://doi.org/10.1101/2020.01.27.921379

% load the preprocessed acceleration, velocity, displacement traces and
% preprocessed EMG
fitsData = load("demoData.mat").fitsData;

% calculate the SRM fits. this will populate the SRMFits table. If you
% would like to calculate the whole thing, which will take a while, unset
% the "EXEMPLARSONLY" flag. Both versions run essentially the same
% functions, but "calculateExemplarSRMFits" has everything for 90°
% perturbations stripped out.

% identify SRM parameters for exemplar patients shown in Figure 2 of the
% manuscript
exemplarPatients = ["bat206" "bat114" "pdf027"];
pertdir = 270;
side = "L";

% isolate only the records described above.
exemplarFitsData = fitsData(ismember(fitsData.patient,exemplarPatients)&ismember(fitsData.pertdir,pertdir)&ismember(fitsData.side,side),:);

SRMFits = calculateExemplarSRMFits(exemplarFitsData);

% plot
plotExemplarSRMFits(SRMFits)

end

function plotExemplarSRMFits(data)
% function plotExemplarSRMFits(data)
%
% plot three example cases

exemplarPatients = ["bat206" "bat114" "pdf027"];
pertdir = 270;
side = "L";

% bat206:
% 22F
% MoCA = 26

% bat114:
% 64M
% MoCA = 26

% pdf027:
% 62F
% pd duration 4.9y
% mds-updrs-iii 55/132
% no freezing

flipTA = false;

XL = [0 1];
YL = [0 1];

% flag to highlight kaPrime in plot
plotComponents = false;

atime = data.atime(1,:);
lookup = atime>=min(XL)&atime<max(XL);
fig = figure;
for pi = 1:length(exemplarPatients)
    ta = data(data.patient==exemplarPatients(pi)&data.pertdir==pertdir&data.side==side&data.mus=="TA",:);
    mg = data(data.patient==exemplarPatients(pi)&data.pertdir==pertdir&data.side==side&data.mus=="MGAS",:);
    
    if flipTA
        s = subplot(2,3,pi);
        xlim(XL)
        ylim(YL)
        
        plot(atime(lookup),mg.e(lookup),'k','linewidth',0.5,'clipping','off');
        plot(atime(lookup),mg.eRecon(lookup),'g','linewidth',0.5,'clipping','off');
        plot(atime(lookup),mg.eRecon(lookup),'k','linewidth',1,'clipping','off');
        
        s = subplot(2,3,pi+3)
        xlim(XL)
        ylim(sort(-1*YL))
        
        plot(atime(lookup),-ta.e(lookup),'k','linewidth',0.5,'clipping','off');
        plot(atime(lookup),-ta.eRecon(lookup),'g','linewidth',0.5,'clipping','off');
        plot(atime(lookup),-ta.ePrimeRecon(lookup),'r','linewidth',0.5,'clipping','off');
        plot(atime(lookup),-ta.eTotalRecon(lookup),'k','linewidth',1,'clipping','off');
    else
        s = subplot(2,3,pi); hold on;
        xlim(XL)
        ylim(YL)
        
        a = area(atime(lookup),mg.eRecon(lookup));
        a.FaceColor = [0 1 0];
        a.EdgeColor = 'none';
        
        plot(atime(lookup),mg.e(lookup),'k','linewidth',0.5,'clipping','off');
        plot(atime(lookup),mg.eRecon(lookup),'k','linewidth',1,'clipping','off');
        
        VAF = rsqr_uncentered(mg.e(lookup)',mg.eRecon(lookup)');
        R2 = rsqr(mg.e(lookup)',mg.eRecon(lookup)');
        VAFstr = "VAF = "+sprintf('%0.2f',VAF);
        R2str = "R2 = "+sprintf('%0.2f',R2);
        txt = text(max(s.XLim),max(s.YLim),[VAFstr;R2str]);
        txt.HorizontalAlignment = 'right';
        txt.VerticalAlignment = 'top';
        
        s = subplot(2,3,pi+3); hold on;
        xlim(XL)
        ylim(YL)
        
        if plotComponents
            componentSub(ta)
        else
            a = area(atime(lookup),ta.eRecon(lookup));
            a.FaceColor = [0 1 0];
            a.EdgeColor = 'none';
            a = area(atime(lookup),ta.ePrimeRecon(lookup));
            a.FaceColor = [1 0 0];
            a.EdgeColor = 'none';
            plot(atime(lookup),ta.e(lookup),'k','linewidth',0.5,'clipping','off');
            plot(atime(lookup),ta.eTotalRecon(lookup),'k','linewidth',1,'clipping','off');
        end
        VAF = rsqr_uncentered(ta.e(lookup)',ta.eTotalRecon(lookup)');
        R2 = rsqr(ta.e(lookup)',ta.eTotalRecon(lookup)');
        VAFstr = "VAF = "+sprintf('%0.2f',VAF);
        R2str = "R2 = "+sprintf('%0.2f',R2);
        txt = text(max(s.XLim),max(s.YLim),[VAFstr;R2str]);
        txt.HorizontalAlignment = 'right';
        txt.VerticalAlignment = 'top';
    end
end

    function componentSub(ta)
        kaPrime = ta.xTotal(5);
        kvPrime = ta.xTotal(6);
        kdPrime = ta.xTotal(7);
        lambdaPrime = ta.xTotal(8);
        
        aPrime = scaleAndDelay(ta.aPrime,lambdaPrime,kaPrime);
        vPrime = scaleAndDelay(ta.vPrime,lambdaPrime,kvPrime);
        dPrime = scaleAndDelay(ta.dPrime,lambdaPrime,kdPrime);
        vdPrime = threshold( kvPrime*channelDelay(ta.vPrime,lambdaPrime,atime) + kdPrime*channelDelay(ta.dPrime,lambdaPrime,atime))
        
        a = area(atime(lookup),ta.eRecon(lookup));
        a.FaceColor = [0 1 0];
        a.EdgeColor = 'none';
        
        a = area(atime(lookup),aPrime(lookup));
        a.FaceColor = [1 0 0];
        a.EdgeColor = 'none';
        
        a = area(atime(lookup),vdPrime(lookup));
        a.FaceColor = [1 0.55 0];
        a.EdgeColor = 'none';
        
        plot(atime(lookup),ta.e(lookup),'k','linewidth',0.5,'clipping','off');
        plot(atime(lookup),ta.eTotalRecon(lookup),'k','linewidth',1,'clipping','off');
        
    end

    function out = scaleAndDelay(in,delay,gain)
        out = threshold(gain*channelDelay(in,delay,atime))
    end

    function out = channelDelay(in,delay,atime)
        out = interp1(atime,in,atime-delay,'linear',0);
    end

    function out = threshold(in)
        out = max(in,0);
    end

end

function SRMFits = calculateExemplarSRMFits(fitsData)

% use a common copy of the atime vector
atime = fitsData.atime(1,:);

removeBackLev = true;
fitDestabMG = true;

% loop. note that randperm here does not affect functionality and is just
% used for debugging purposes (in order to prevent having to go through all
% of the TA records prior to doing MG.)
for idx = randperm(size(fitsData,1))
    e = fitsData.e(idx,:);
    
    if removeBackLev
        backLev = nanmean(e(atime<0.05));
        e = e - backLev;
    end
    if fitsData.mus(idx)=="TA"&fitsData.dir(idx)=="B"
        
        % identify braking response in TA in backward perturbations
        predictorsBraking = [fitsData.a(idx,:); fitsData.v(idx,:); fitsData.d(idx,:)];
        [fitsData.x(idx,[1:4]), fitsData.eRecon(idx,:) fitsData.fit(idx,:)] = fitSRM(e,atime,predictorsBraking);
        
        % identify destabilizing response in TA in backward perturbations
        predictorsDestabilizing = [fitsData.aPrime(idx,:); fitsData.vPrime(idx,:); fitsData.dPrime(idx,:)];
        [fitsData.xPrime(idx,[5:8]), fitsData.ePrimeRecon(idx,:) fitsData.fitPrime(idx,:)] = fitDestabilizingSRM(e,atime,predictorsDestabilizing);
        
        % combine them into the initial guess for the final optimization
        X0Total = [fitsData.x(idx,[1:4]) fitsData.xPrime(idx,[5:8])];
        predictorsTotal = [predictorsBraking; predictorsDestabilizing];
        
        [fitsData.xTotal(idx,:), fitsData.eTotalRecon(idx,:) fitsData.fitTotal(idx,:)] = fitTotalSRM(e,atime,predictorsTotal,X0Total);
        if removeBackLev
            fitsData.eTotalRecon(idx,:) = fitsData.eTotalRecon(idx,:) + backLev;
        end
    elseif fitsData.mus(idx)=="MGAS"&fitsData.dir(idx)=="B"
        % identify braking response in MGAS in backward perturbations
        predictorsBraking = [fitsData.a(idx,:); fitsData.v(idx,:); fitsData.d(idx,:)];
        [fitsData.x(idx,[1:4]), fitsData.eRecon(idx,:) fitsData.fit(idx,:)] = fitSRM(e,atime,predictorsBraking);
        if removeBackLev
            fitsData.eRecon(idx,:) = fitsData.eRecon(idx,:) + backLev;
        end
    end
end

SRMFits = fitsData;
end

function SRMFits = calculateSRMFits(fitsData)

% use a common copy of the atime vector
atime = fitsData.atime(1,:);

removeBackLev = true;
fitDestabMG = true;

% loop. note that randperm here does not affect functionality and is just
% used for debugging purposes (in order to prevent having to go through all
% of the TA records prior to doing MG.)
for idx = randperm(size(fitsData,1))
    e = fitsData.e(idx,:);
    
    if removeBackLev
        backLev = nanmean(e(atime<0.05));
        e = e - backLev;
    end
    if fitsData.mus(idx)=="TA"&fitsData.dir(idx)=="B"
        
        % identify braking response in TA in backward perturbations
        predictorsBraking = [fitsData.a(idx,:); fitsData.v(idx,:); fitsData.d(idx,:)];
        [fitsData.x(idx,[1:4]), fitsData.eRecon(idx,:) fitsData.fit(idx,:)] = fitSRM(e,atime,predictorsBraking);
        
        % identify destabilizing response in TA in backward perturbations
        predictorsDestabilizing = [fitsData.aPrime(idx,:); fitsData.vPrime(idx,:); fitsData.dPrime(idx,:)];
        [fitsData.xPrime(idx,[5:8]), fitsData.ePrimeRecon(idx,:) fitsData.fitPrime(idx,:)] = fitDestabilizingSRM(e,atime,predictorsDestabilizing);
        
        % combine them into the initial guess for the final optimization
        X0Total = [fitsData.x(idx,[1:4]) fitsData.xPrime(idx,[5:8])];
        predictorsTotal = [predictorsBraking; predictorsDestabilizing];
        
        [fitsData.xTotal(idx,:), fitsData.eTotalRecon(idx,:) fitsData.fitTotal(idx,:)] = fitTotalSRM(e,atime,predictorsTotal,X0Total);
        if removeBackLev
            fitsData.eTotalRecon(idx,:) = fitsData.eTotalRecon(idx,:) + backLev;
        end
    elseif fitsData.mus(idx)=="TA"&fitsData.dir(idx)=="F"
        % identify braking response in TA in forward perturbations
        predictorsBraking = [fitsData.a(idx,:); fitsData.v(idx,:); fitsData.d(idx,:)];
        [fitsData.x(idx,[1:4]), fitsData.eRecon(idx,:) fitsData.fit(idx,:)] = fitSRM(e,atime,predictorsBraking);
        if removeBackLev
            fitsData.eRecon(idx,:) = fitsData.eRecon(idx,:) + backLev;
        end
    elseif fitsData.mus(idx)=="MGAS"&fitsData.dir(idx)=="B"
        % identify braking response in MGAS in backward perturbations
        predictorsBraking = [fitsData.a(idx,:); fitsData.v(idx,:); fitsData.d(idx,:)];
        [fitsData.x(idx,[1:4]), fitsData.eRecon(idx,:) fitsData.fit(idx,:)] = fitSRM(e,atime,predictorsBraking);
        if removeBackLev
            fitsData.eRecon(idx,:) = fitsData.eRecon(idx,:) + backLev;
        end
    elseif fitsData.mus(idx)=="MGAS"&fitsData.dir(idx)=="F"
        % identify braking response in MGAS in forward perturbations
        predictorsBraking = [fitsData.a(idx,:); fitsData.v(idx,:); fitsData.d(idx,:)];
        [fitsData.x(idx,[1:4]), fitsData.eRecon(idx,:) fitsData.fit(idx,:)] = fitSRM(e,atime,predictorsBraking);
        
        if fitDestabMG
            % identify destabilizing response in MG in forward perturbations
            predictorsDestabilizing = [fitsData.aPrime(idx,:); fitsData.vPrime(idx,:); fitsData.dPrime(idx,:)];
            [fitsData.xPrime(idx,[5:8]), fitsData.ePrimeRecon(idx,:) fitsData.fitPrime(idx,:)] = fitDestabilizingSRM(e,atime,predictorsDestabilizing);
            
            % combine them into the initial guess for the final optimization
            X0Total = [fitsData.x(idx,[1:4]) fitsData.xPrime(idx,[5:8])];
            predictorsTotal = [predictorsBraking; predictorsDestabilizing];
            
            [fitsData.xTotal(idx,:), fitsData.eTotalRecon(idx,:) fitsData.fitTotal(idx,:)] = fitTotalSRM(e,atime,predictorsTotal,X0Total);
            if removeBackLev
                fitsData.eTotalRecon(idx,:) = fitsData.eTotalRecon(idx,:) + backLev;
            end
        end
        if removeBackLev
            fitsData.eRecon(idx,:) = fitsData.eRecon(idx,:) + backLev;
        end
    end
end

SRMFits = fitsData;

end

function out = assembleChannel(signals,gains,delay,atime)
out = channelDelay(threshold(signals'*gains'),delay,atime);
end

function out = assembleChannelComponents(signals,gains,delay,atime,THRESH)

% note that we cannot just matrix multiply because we have to delay each signal
out = [];
if THRESH
    for i = 1:length(gains)
        out(i,:) = channelDelay(threshold(signals(i,:)*gains(i)),delay,atime);
    end
else
    for i = 1:length(gains)
        out(i,:) = channelDelay(signals(i,:)*gains(i),delay,atime);
    end
end
end

function out = assembleTwoChannels(signals1,gains1,delay1,signals2,gains2,delay2,atime)

out1 = channelDelay(threshold(signals1'*gains1'),delay1,atime);
out2 = channelDelay(threshold(signals2'*gains2'),delay2,atime);

out = threshold(out1+out2);

    function out = channelDelay(in,delay,atime)
        out = interp1(atime,in,atime-delay,'linear',0);
    end

    function out = threshold(in)
        out = max(in,0);
    end

end

function out = channelDelay(in,delay,atime)
out = interp1(atime,in,atime-delay,'linear',0);
end

function [xBraking, eBraking, fitsBraking] = fitBrakingSRM(e,atime,predictors)
% fit "braking response" at end of perturbation.

% load optimization options. these include overall cost function options as
% well as direct options for the Matlab optimizer.
optimizationParameters = loadOptimizationParameters();

options = optimoptions('fmincon',...
    'display',optimizationParameters.display,...
    'TolX',optimizationParameters.TolX,...
    'TolFun',optimizationParameters.TolFun,...
    'MaxFunEvals',optimizationParameters.MaxFunEvals...
    );

% 60-250 ms search range was welch 2008
X0 = [10.0 0.01 0.01 0.150];
UB = [15.0 0.04 0.04 0.250];
% UB = [15.0 0.04 0.04 0.200];
% UB = [15.0 0.04 0.04 0.175];
% UB = [15.0 0.02 0.02 0.175];
LB = [ 0.0 0.00 0.00 0.060];
% LB = [ 0.0 0.00 0.00 0.050];
% LB = [ 0.0 0.00 0.00 0.080];
% LB = [ 0.0 0.00 0.00 0.100];
gainFlag = [true true true false];

fixBurstGain = true;
if fixBurstGain
    burstGain = max(e(1,atime>0.650&atime<0.800))/max(predictors(1,:));
    [X0(1),UB(1)] = deal(max(min(burstGain,15),0));
    LB(1) = 0.8 * X0(1);
end

[X,FVAL,EXITFLAG] = fmincon(@(X) jigsawPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters),X0,[],[],[],[],LB,UB,[],options);

xBraking = X;

eBraking = assembleChannel(predictors,X(gainFlag),X(~gainFlag),atime);
fitParameters = modelIPI('X',X,'eRecon',eBraking,'e',e,'optimizationParameters',optimizationParameters,'gainWeight',double(gainFlag));

fitsBraking(1) = fitParameters.r2;
fitsBraking(2) = fitParameters.r2u;

fitsBraking(isnan(fitsBraking)) = 0;

end

function [xDestabilizing, eDestabilizing, fitsDestabilizing] = fitDestabilizingSRM(e,atime,predictors)
% fit "Destabilizing response" at end of perturbation.

% load optimization options. these include overall cost function options as
% well as direct options for the Matlab optimizer.
optimizationParameters = loadOptimizationParameters();

options = optimoptions('fmincon',...
    'display',optimizationParameters.display,...
    'TolX',optimizationParameters.TolX,...
    'TolFun',optimizationParameters.TolFun,...
    'MaxFunEvals',optimizationParameters.MaxFunEvals...
    );

% note slightly different gain limits for destabilizing SRM TA gain
X0 = [10.0 0.01 0.01 0.140];
% X0 = [10.0 0.01 0.01 0.150];
UB = [15.0 0.04 0.04 0.210];
% UB = [15.0 0.04 0.04 0.250];
% UB = [15.0 0.04 0.04 0.200];
% UB = [15.0 0.04 0.04 0.175];
% UB = [15.0 0.02 0.02 0.175];
LB = [ 0.0 0.00 0.00 0.090];
% LB = [ 0.0 0.00 0.00 0.060];
% LB = [ 0.0 0.00 0.00 0.050];
% LB = [ 0.0 0.00 0.00 0.080];
% LB = [ 0.0 0.00 0.00 0.100];
gainFlag = [true true true false];

fixBurstGain = true;
if fixBurstGain
    burstGain = max(e(1,atime>0.150&atime<0.275))/max(predictors(1,:));
    % burstGain = max(e(1,atime>0.150&atime<0.300))/max(predictors(1,:));
    [X0(1),UB(1)] = deal(max(min(burstGain,15),0));
    LB(1) = 0.9 * X0(1);
end

[X,FVAL,EXITFLAG] = fmincon(@(X) jigsawPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters),X0,[],[],[],[],LB,UB,[],options);

xDestabilizing = X;

eDestabilizing = assembleChannel(predictors,X(gainFlag),X(~gainFlag),atime);
fitParameters = modelIPI('X',X,'eRecon',eDestabilizing,'e',e,'optimizationParameters',optimizationParameters,'gainWeight',double(gainFlag));

fitsDestabilizing(1) = fitParameters.r2;
fitsDestabilizing(2) = fitParameters.r2u;
fitsDestabilizing(isnan(fitsDestabilizing)) = 0;

end

function [xBraking, eBraking, fitsBraking] = fitSRM(e,atime,predictors)
% fit "braking response" at end of perturbation.

% load optimization options. these include overall cost function options as
% well as direct options for the Matlab optimizer.
optimizationParameters = loadOptimizationParameters();

options = optimoptions('fmincon',...
    'display',optimizationParameters.display,...
    'TolX',optimizationParameters.TolX,...
    'TolFun',optimizationParameters.TolFun,...
    'MaxFunEvals',optimizationParameters.MaxFunEvals...
    );

% 60-250 ms search range as in welch 2008
X0 = [10.0 0.01 0.01 0.150];
UB = [15.0 0.04 0.04 0.250];
% UB = [15.0 0.04 0.04 0.200];
% UB = [15.0 0.04 0.04 0.175];
% UB = [15.0 0.02 0.02 0.175];
LB = [ 0.0 0.00 0.00 0.060];
% LB = [ 0.0 0.00 0.00 0.050];
% LB = [ 0.0 0.00 0.00 0.080];
% LB = [ 0.0 0.00 0.00 0.100];
gainFlag = [true true true false];

fixBurstGain = true;
if fixBurstGain
    burstGain = max(e(1,atime>0.150&atime<0.800))/max(predictors(1,:));
    [X0(1),UB(1)] = deal(max(min(burstGain,15),0));
    LB(1) = 0.9 * X0(1);
end

[X,FVAL,EXITFLAG] = fmincon(@(X) jigsawPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters),X0,[],[],[],[],LB,UB,[],options);

xBraking = X;

eBraking = assembleChannel(predictors,X(gainFlag),X(~gainFlag),atime);
fitParameters = modelIPI('X',X,'eRecon',eBraking,'e',e,'optimizationParameters',optimizationParameters,'gainWeight',double(gainFlag));

fitsBraking(1) = fitParameters.r2;
fitsBraking(2) = fitParameters.r2u;
fitsBraking(isnan(fitsBraking)) = 0;

end

function [xTotal, eTotal,fitsTotal] = fitTotalSRM(e,atime,predictors,X0)

% load optimization options. these include overall cost function options as
% well as direct options for the Matlab optimizer.
optimizationParameters = loadOptimizationParameters();

options = optimoptions('fmincon',...
    'display',optimizationParameters.display,...
    'TolX',optimizationParameters.TolX,...
    'TolFun',optimizationParameters.TolFun,...
    'MaxFunEvals',optimizationParameters.MaxFunEvals...
    );

% allow the gains to vary within +/-20%, allow the delays to vary within +/- 20 msec. note that lower bound for gains must be positive.
gainFlag = [true true true false true true true false];
UB(gainFlag) = 1.1*X0(gainFlag);
LB(gainFlag) = max(0.9*X0(gainFlag),0);
[UB(~gainFlag),LB(~gainFlag)] = deal(X0(~gainFlag));

% allow the final delays to vary within +/- 10sec
[UB(~gainFlag),LB(~gainFlag)] = deal(X0(~gainFlag));
UB(~gainFlag) = min(X0(~gainFlag)+0.010,0.250);
LB(~gainFlag) = max(X0(~gainFlag)-0.010,0.060);

% add a very small offset to improve convergence when LB and UB are very
% close.
UB((UB-LB)<1e-6) = UB((UB-LB)<1e-6)+1e-6;

[X,FVAL,EXITFLAG] = fmincon(@(X) jigsawTwoChannelPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters),X0,[],[],[],[],LB,UB,[],options);

xTotal = X;
eTotal = assembleTwoChannels(predictors(1:3,:),X(1:3),X(4),predictors(4:6,:),X(5:7),X(8),atime);

fitsTotal(1) = rsqr(e',eTotal');
fitsTotal(2) = rsqr_uncentered(e',eTotal');

end

function [xBraking, eBraking, fitsBraking] = fitTraditionalSRM(e,atime,predictors)
% fit "braking response" at end of perturbation.

% load optimization options. these include overall cost function options as
% well as direct options for the Matlab optimizer.
optimizationParameters = loadOptimizationParameters();

options = optimoptions('fmincon',...
    'display',optimizationParameters.display,...
    'TolX',optimizationParameters.TolX,...
    'TolFun',optimizationParameters.TolFun,...
    'MaxFunEvals',optimizationParameters.MaxFunEvals...
    );

% 60-250 ms search range as in welch 2008
X0 = [10.0 0.01 0.01 0.150];
UB = [15.0 0.04 0.04 0.250];
% UB = [15.0 0.04 0.04 0.200];
% UB = [15.0 0.04 0.04 0.175];
% UB = [15.0 0.02 0.02 0.175];
LB = [ 0.0 0.00 0.00 0.060];
% LB = [ 0.0 0.00 0.00 0.050];
% LB = [ 0.0 0.00 0.00 0.080];
% LB = [ 0.0 0.00 0.00 0.100];
gainFlag = [true true true false];

fixBurstGain = true;
if fixBurstGain
    burstGain = max(e(1,atime>0.150&atime<0.300))/max(predictors(1,:));
    [X0(1),UB(1)] = deal(max(min(burstGain,15),0));
    LB(1) = 0.9 * X0(1);
end

[X,FVAL,EXITFLAG] = fmincon(@(X) jigsawPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters),X0,[],[],[],[],LB,UB,[],options);

xBraking = X;

eBraking = assembleChannel(predictors,X(gainFlag),X(~gainFlag),atime);
fitParameters = modelIPI('X',X,'eRecon',eBraking,'e',e,'optimizationParameters',optimizationParameters,'gainWeight',double(gainFlag));

fitsBraking(1) = fitParameters.r2;
fitsBraking(2) = fitParameters.r2u;
fitsBraking(isnan(fitsBraking)) = 0;

end

function PI = jigsawPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters)

% calculate reconstruction
eRecon = assembleChannel(predictors,X(gainFlag),X(~gainFlag),atime);

% calculate performance index
fitParameters = modelIPI('X',X,'eRecon',eRecon,'e',e,'optimizationParameters',optimizationParameters,'gainWeight',double(gainFlag));

PI = fitParameters.PI;

end

function PI = jigsawTwoChannelPassthrough(X,predictors,gainFlag,atime,e,optimizationParameters)

% calculate reconstruction
eRecon = assembleTwoChannels(predictors(1:3,:),X(1:3),X(4),predictors(4:6,:),X(5:7),X(8),atime);

% calculate performance index
fitParameters = modelIPI('X',X,'eRecon',eRecon,'e',e,'optimizationParameters',optimizationParameters,'gainWeight',double(gainFlag));

PI = fitParameters.PI;
end

function out = loadOptimizationParameters()
J_SE = 1;
J_ME = 1;
J_GM = 1;
display = "iter";
TolX = 1e-9;
MaxFunEvals = 1e+5;
TolFun = 1e-7;
out = table(J_SE,J_ME,J_GM,display,TolX,MaxFunEvals,TolFun);
end

function fitParameters = modelIPI(varargin)

p = inputParser;
p.addOptional('X',[]);
p.addOptional('e',[]);
p.addOptional('eRecon',[]);
p.addOptional('gainWeight',[]);
p.addOptional('optimizationParameters',[]);
p.parse(varargin{:});

X = p.Results.X;
e = p.Results.e;
e_recon = p.Results.eRecon;
optimizationParameters = p.Results.optimizationParameters;
gainWeight = p.Results.gainWeight;

J_SE = optimizationParameters.J_SE;
J_ME = optimizationParameters.J_ME;
J_GM = optimizationParameters.J_GM;

fitParameters.r2 = rsqr(e',e_recon');
fitParameters.r2u = rsqr_uncentered(e',e_recon');

% Compute the errors between recorded and simulated EMG signals. We are
% minimizing a function of these error values:
% error   = recorded  - simulated
e_error = e - e_recon;
e_error(isnan(e_error)) = [];

% ### J = squared error + "min-max error" = terminal cost
% ### J = L1 + L2:
% L1 = squared error. L1 composes the bulk of the terminal cost. It is a
% linear combination of the squares of the error terms.
fitParameters.L1 = J_SE*(e_error*e_error');

% L2 = max of the abs of the error terms. L2 penalizes large deviations.
% Little terminal cost is contributed by L2, but L2 appears to "smooth out"
% the J manifold and make the optimization converge more consistently.
fitParameters.L2 = J_ME*max(abs(e_error));

% L3 = magnitude of the gain values. this term is just designed to zero out
% noncontributing gains, and improve convergence.
fitParameters.L3 = J_GM*(X.*gainWeight)*(X.*gainWeight)';

% PI is the sum of these individual costs. The function fmincon attempts to
% minimize J.
fitParameters.PI = fitParameters.L1 + fitParameters.L2 + fitParameters.L3;

end

function ursqr = rsqr_uncentered(data,data_rec)
% This function calculates the uncetered correlation coefficient using "Cluster" method.
%
% Syntax:   r_sqr = rsqr_uncentered(data,data_rec)
%
% Input:
% data      Array   matrix of observed data  (e.g., data = [observations x channels])
% data_rec  Array   matrix of reconstructed/predicted data (e.g., data_rec = [observations x channels])
%
% Output:
% ursqr     Array   row vector of uncentered correlation coefficients
%
% Created: May 24, 2006 (Gelsy Torres-Oviedo)
% Last Modified: July 10, 2006 (Torrence Welch)
% Last Modification: fix ursqr calculation
% Last Modified: June 7, 2016 (J. Lucas McKay)
% Last Modification:
%   remove transpose to expect data in columns
%   remove empty values

% Zar book method p. 334
nvars = size(data,2);
ursqr = nan(1,nvars);
for i = 1:nvars
    X = [data(:,i) data_rec(:,i)];
    X(isnan(X(:,1))|isnan(X(:,2)),:) = [];
    ursqr(i) = sum(prod(X,2))^2 / (sum(X(:,1).^2)*sum(X(:,2).^2)); %regression sum of squares/total sum of squares
end

end

function rsqr = rsqr(data,data_rec)
% calculates coefficient of determination (R-squared) according to
% methodology suggested by the Mathworks
%
% https://www.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
% accessed 2017 07 07
%
% Created: June 7, 2017 (J. Lucas McKay)
% last edit: 3/20/19: corrected an error in indexing the STATS structure
% that failed if the number of data columns was >1.

[nrows,ncols] = size(data);
rsqr = nan(1,ncols);

% method = categorical(cellstr('fitlm'));
method = categorical(cellstr('regress'));

if method == 'fitlm'
    for i = 1:ncols
        X = data_rec(:,i);
        y = data(:,i);
        mdl = fitlm(X,y);
        rsqr(i) = mdl.Rsquared.Ordinary;
    end
elseif method == 'regress'
    for i = 1:ncols
        X = [data_rec(:,i) ones(size(data_rec(:,i)))];
        y = data(:,i);
        [B,~,~,~,stats] = regress(y,X);
        rsqr(i) = stats(1);
    end
end

% the following instances demonstrate that the regress method is
% approximatley 15x faster than the fitlm method.
if false
    f = [];
    r = [];
    for j = 1:1000
        tic
        X = data_rec(:,i);
        y = data(:,i);
        mdl = fitlm(X,y);
        rsqr(i) = mdl.Rsquared.Ordinary;
        f(end+1) = toc;
        
        tic
        X = [data_rec(:,i) ones(size(data_rec(:,i)))];
        y = data(:,i);
        [B,~,~,~,stats] = regress(y,X);
        rsqr(i) = stats(1);
        r(end+1) = toc;
    end
    mean(f)/mean(r)
end

end

function a_out = stictionTemplate(varargin)

p = inputParser;
p.addOptional('a',[]);

% p.addOptional('t',[0:(1/1080):(0.8-(1/1080))]);
p.addOptional('t',[-0.5:(9.2593e-04):0.9991]);
p.addOptional('ts1',0.000);
p.addOptional('ts2',0.075); % initial burst; search range 50-100 ms
p.addOptional('ts3',0.450);
p.addOptional('ts4',0.450);
% p.addOptional('ts5',0.525); % deceleration burst; search range 500-550 ms
p.addOptional('ts5',0.650); % deceleration burst; search range 500-550 ms
p.addOptional('tau',0.005);
p.parse(varargin{:});

a = p.Results.a;
t = p.Results.t;

ts1 = p.Results.ts1;
ts2 = p.Results.ts2;
ts3 = p.Results.ts3;
ts4 = p.Results.ts4;
ts5 = p.Results.ts5;
tau = p.Results.tau;

a_out = nan(size(a));

% identify indices here
bkgdInds = t<ts1;
burst1Inds = findInds(ts1,ts2);
silence1Inds = findInds(ts2,ts3);
burst2Inds = findInds(ts4,ts5);
silence2Inds = find(t>ts5);

a_out(:,bkgdInds) = 0;

for i = 1:size(a_out,1)
    a_orig = a(i,:);
    a_new = nan(size(a_orig));
    a_new(bkgdInds) = 0;
    
    a_new(burst1Inds) = a_orig(burst1Inds) - a_orig(burst1Inds(1));
    silence1Time = t(silence1Inds) - t(silence1Inds(1));
    a_new(silence1Inds) = a_new(burst1Inds(end))*exp(-(1/tau)*silence1Time);
    
    a_new(burst2Inds) = a_orig(burst2Inds) - a_orig(burst2Inds(1));
    silence2Time = t(silence2Inds) - t(silence2Inds(1));
    a_new(silence2Inds) = a_new(burst2Inds(end))*exp(-(1/tau)*silence2Time);
    
    a_out(i,:) = a_new;
end

    function out = findInds(tstart,tstop)
        out = [time2ind(t,tstart):time2ind(t,tstop)];
    end

end

function out = threshold(in)
out = max(in,0);
end








