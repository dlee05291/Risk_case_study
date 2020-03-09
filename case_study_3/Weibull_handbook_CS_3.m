% Variables Info.(alphabetical order) -------------------------------------
% Beta: Weibull shape parameter
% DeltaFcount: Failure counts during each month
% DeltaNew: Number of renewed items(= FiCount)
% DeltaPOF: Probability of failure for each month
% Eta: Weibull location parameter
% FF: Failure forecast function
% FfCount: Total failure numbers(initial + renewed)
% FFtime: Operation time for whole months
% FFtimeUnit: Operation time after certain month
% FiCount: Failure numbers of initial items
% FidxPerMonth: Logical indices of failed unit for MCS
% idxnp: Logical indices for resampling of non-physical samples for MCS
% MCSdeltaNEW: Number of renewed items
% MCSFfCount: Total failure numbers(initial + renewed)
% MCSFnCount: Number of failured items among renewed items
% MCSiLife: Life distribution of initial items
% MCSnLife: Life distribution of renewed items
% MCSSfCount: System failure counts
% Month: Months to forecast(constant)
% POFi: Failure probability of initial items
% POFn: Failure probability of renewed items
% Samples: Monte Carlo simulation sample number(constant)
% SfCount: System total failure numbers
% Tsn: Time since new for each items(same as population number)
% TTF: Weibully distributed life generation function
% Usage: Operation time per month

% Execute -----------------------------------------------------------------
clear all; close all; clc;

% Generate Items
TempA = ones(10,1) * 10;
TempB = ones(10,1) * 20;
Tsn = cat(1,TempA,TempB);
clear TempA TempB;

% Input variables
Usage = 25;                             % operation time per month
Month = 60;                             % months to forecast
Samples = 1e5;                          % sample number for MCS

% Constants & Outputs from input variables
Beta = [0.5 1 2 6];
Eta = [30000 1200 700 600];
FFtime = (1:Month) * Usage;             % Failure Forecasting time

% Failure forcast function
FF = @(t0, t1, s) (wblcdf(t1, Eta(s), Beta(s)) - ...
    wblcdf(t0, Eta(s), Beta(s)))./ ...
    (1 - wblcdf(t0, Eta(s), Beta(s)));
% s: failure mode

% Weibully distributed life generation function
TTF = @(FM, SN) wblrnd(Eta(FM), Beta(FM), [SN, 1]);
% FM: failure mode, SN: sample number

% Failure forecast with renewal -------------------------------------------

% Initialize variables regarding failure items
FiCount = zeros(numel(FFtime), 1);          % failures from initial items
DeltaNew = zeros(numel(FFtime), 1);         % renewed items(= FiCount)
DeltaFcount = zeros(numel(FFtime), 1);      % failures from DeltaNew

% Failure forecast based on failure forecast function ---------------------
for j = 1:numel(Beta)                       % j means no. of failure modes
    
    for i = 1:Month
        % Number of failure for initial unit(month)
        FiCount(i,j) = sum(FF(Tsn, Tsn + FFtime(i), j));
        
        % Number of failure for renewed unit
        if i == 1
            % Number of new units introduced
            DeltaNew(i,j) = FiCount(i,j);
            
        else
            % Number of new units introduced
            DeltaNew(i,j) = FiCount(i,j) - FiCount(i-1,j);
            
            % Increase in probability of failure for current step
            DeltaPOF = FF(0, flip(FFtime(1:i-1)), j);
        
            % Increase in number of failure for current step.
            DeltaFcount(i,j) = sum(DeltaPOF(:).*DeltaNew(1:i-1,j));
        end
    end
end

% Get failure counts (initial + renewed)
FfCount = FiCount + DeltaFcount;

% Get system failure counts
SfCount = sum(FfCount,2);

% Failure forecast with MCS -----------------------------------------------

POFi = zeros(numel(Beta),numel(Tsn), numel(FFtime));  % POF of initial items
MCSdeltaNEW = zeros(numel(Beta),numel(Tsn), numel(FFtime)); 

for k = 1:numel(Beta)                           % k: no. of failure modes
    
    % Initialize variables for POF and no. of renewd item    
    
    for j = 1:numel(Tsn)                        % j; no. of items
        % Generate Weibull initial life distribution for each failure mode
        MCSiLife(:,k) = TTF(k, Samples);
        
        % Resample for non-physical samples(less than TSN)
        while sum(MCSiLife < Tsn(j)) ~= 0
            idxnp = MCSiLife < Tsn(j);
            MCSiLife(idxnp) = TTF(k, sum(idxnp));
        end
        
            for i = 1:Month
                % Failure forecast for initial items-----------------------
                % Forecast time after i month. Constant for comparison.
                FFtimeUnit = Tsn(j) + FFtime(i);
                
                % Indices for failed unit.
                FidxPerMonth = MCSiLife(:,k) < FFtimeUnit;
                
                % Get POF for initial unit.
                POFi(k,j,i) = sum(FidxPerMonth)/Samples;
                
                % Failure forecast for renewal items ----------------------
                if i == 1
                    % Set renewal item POF for first month as zero.
                    POFn(k,j,i) = 0;
                    
                    % Generate renewal life distribution.
                    MCSnLife(:,k) = TTF(k, Samples);
                    
                else 
                    % Calculate failure rate during one month(i-1 to i)
                    % No. of renewal items during one month
                    MCSdeltaNEW(k,j,i) = (POFi(k,j,i) - POFi(k,j,i-1));
                    
                    % Calculate POF during one month
                    FidxPerMonth = MCSnLife(:,k) < FFtime(i-1);
                    POFn(k,j,i) = (sum(FidxPerMonth)/Samples);
                    
                    % Failure no. of renewed items
                    MCSFnCount(k,j,i) = sum(flip(squeeze(POFn(k,j,1:i))) .* squeeze(MCSdeltaNEW(k,j,1:i)));
                    
                end
            end
    end
end

% Get failure counts (initial + renewed)
MCSFfCount = POFi + MCSFnCount;

% Get system failure counts
MCSSfCount = sum(MCSFfCount,2);

% Plot --------------------------------------------------------------------

% Monthly failure rhythm
figure(2); hold on;
plot([FfCount(1,1); FfCount(2:end,1) - FfCount(1:end-1,1)], 's');
plot([FfCount(1,2); FfCount(2:end,2) - FfCount(1:end-1,2)], 'd');
plot([FfCount(1,3); FfCount(2:end,3) - FfCount(1:end-1,3)], 'o');
plot([FfCount(1,4); FfCount(2:end,4) - FfCount(1:end-1,4)], '*');
plot([SfCount(1); SfCount(2:end) - SfCount(1:end-1)], '+');
grid on
title('Monthly Failure Rhythm');
xlabel('Calandar Time (months)'); ylabel('Failures per Month');
legend({'Control', 'FOD', 'Bearing', 'Turbine', 'System Failure'},'Location','best');

% Cumulative failures
figure(1); hold on;
plot(FfCount(:,1));
plot(FfCount(:,2));
plot(FfCount(:,3));
plot(FfCount(:,4));
plot(SfCount);

plot(permute(sum(squeeze(MCSFfCount(1,:,:)),1),[2 1]), 's');
plot(permute(sum(squeeze(MCSFfCount(2,:,:)),1),[2 1]), 'd');
plot(permute(sum(squeeze(MCSFfCount(3,:,:)),1),[2 1]), 'o');
plot(permute(sum(squeeze(MCSFfCount(4,:,:)),1),[2 1]), '*');
plot(sum(permute(squeeze(MCSSfCount), [2 1]),2), '+');
grid on

axis([0, 60, 0, 140]);
title('Cumulative System Failures');
xlabel('Calandar Time (months)'); ylabel('Cumulative Failures');
legend({'Control', 'FOD', 'Bearing', 'Turbine', 'System Failure', ...
    'MCS-Ctrl', 'MCS-FOD', 'MCS-Bear', 'MCS-Turb', 'MCS-System'},'Location','best');
