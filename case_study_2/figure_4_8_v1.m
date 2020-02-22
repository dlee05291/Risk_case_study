% Figure 4-8. Bleed System With Renewal
%
% Revision history
% 022020 LDY Code was created. Renewal was accomplished by adding new units
%            to suspension fleet when POF of individual unit exceed control
%            level.

clear; clc;

% Constant
usageRate = 50;    % Operation time per month, [h/m]
nsample = 1e6;     % Number of samples for MCS

% Variable
pofCtrl = 0.5;    % POF criteria for introducing new unit

% Months and times for failure forcast
mff = 1:60;
tff = mff*usageRate;

% Load data.
data = csvread('case_study_2_data_for_d.csv');
tsn = data(:, 1);
idxFailed = data(:, 2);

% Replace failed part w/ new one.
%tsn(logical(idxFailed)) = 0;

% Consider only suspension data.
tsn(logical(idxFailed)) = [];

% Number of units
nunit = numel(tsn);

% Parameters
etaParam = 2004.45;
betaParam = 5.239307;

% Without renewal (conventional) ------------------------------------------
% Model for failure forcast
FF = @(t0, t1) (wblcdf(t1, etaParam, betaParam) - ...
    wblcdf(t0, etaParam, betaParam))./ ...
    (1 - wblcdf(t0, etaParam, betaParam));

% Initialize variables for the number of failed unit.
nfailed = zeros(numel(mff), 1);

% Do calculation.
for i = 1:numel(mff)
    nfailed(i) = sum(FF(tsn, tsn + tff(i)));
end

% Get failure forecast.
ffc1 = nfailed;

% With renewal (conventional) ---------------------------------------------
% Initialize variables for the number of failed unit.
nfailed1 = zeros(nunit, numel(tff));
nfailed2 = zeros(nunit, numel(tff));

% Do calculation.
for i = 1:nunit
    % Default renewal flag
    flagRenewal = 0;
    
    for j = 1:numel(tff)
        % Number of failure for initial unit
        nfailed1(i, j) = FF(tsn(i), tsn(i) + tff(j));
        
        % Number of failure for renewed unit
        if flagRenewal == 1
            nfailed2(i, j) = FF(0, tff(j) - tff(idxRenewal));
        end
        
        if flagRenewal == 0 && nfailed1(i, j) > pofCtrl
            % Update renewal flag.
            flagRenewal = 1;
            
            % Save index at renewal.
            idxRenewal = j;
        end
    end
end

% Get failure forecast.
ffc2 = sum(nfailed1 + nfailed2, 1)';

% Without renewal (MCS) ---------------------------------------------------
% Initialize variable for POF.
pof = zeros(nunit, numel(tff));

% Time to failure model
TTF = @(N) wblrnd(etaParam, betaParam, [N, 1]);
%TTF = @(N) etaParam*(log(1./(1 - rand(N, 1)))).^(1/betaParam);

% Do MCS.
for i = 1:nunit
    % Time to failure
    ttf = TTF(nsample);
    
    % Resample for non-physical samples.
    while sum(ttf < tsn(i)) ~= 0
        idxnp = ttf < tsn(i);
        ttf(idxnp) = TTF(sum(idxnp));
    end
    
    for j = 1:numel(tff)
        % Current time for failure forecast
        tffUnit = tsn(i) + tff(j);
        
        % Indices for failed unit
        idxFailedPerUnit = ttf < tffUnit;
        
        % Get the number of failures.
        pof(i, j) = sum(idxFailedPerUnit)/nsample;
    end
end

% Get failure forecast.
ffmc1 = sum(pof, 1)';

% With renewal (MCS) ---------------------------------------------------
% Initialize variable for POF.
pof1 = zeros(nunit, numel(tff));
pof2 = zeros(nunit, numel(tff));

% Time to failure model
TTF = @(N) wblrnd(etaParam, betaParam, [N, 1]);

% Do MCS.
for i = 1:nunit
    % Time to failure
    ttf1 = TTF(nsample);
    
    % Resample for non-physical samples.
    while sum(ttf1 < tsn(i)) ~= 0
        idxnp = ttf1 < tsn(i);
        ttf1(idxnp) = TTF(sum(idxnp));
    end
    
    % Default renewal flag
    flagRenewal = 0;
    
    for j = 1:numel(tff)        
        % Current time for failure forecast
        tffUnit1 = tsn(i) + tff(j);
        
        % Indices for failed unit
        idxFailedPerUnit = ttf1 < tffUnit1;
        
        % Get POF for initial unit.
        pof1(i, j) = sum(idxFailedPerUnit)/nsample;
        
        % Get POF for renewed unit
        if flagRenewal == 1
            tffUnit2 = tff(j) - tff(idxRenewal);
            idxFailedPerUnit2 = ttf2 < tffUnit2;
            pof2(i, j) = sum(idxFailedPerUnit2)/nsample;
        end
        
        % Renewal
        if flagRenewal == 0 && pof1(i, j) > pofCtrl
            % Time to failure
            ttf2 = TTF(nsample);
            
            % Update renewal flag.
            flagRenewal = 1;
            
            % Save index at renewal.
            idxRenewal = j;
        end
    end
end

% Get failure forecast.
ffmc2 = sum(pof1 + pof2, 1)';

% Plot --------------------------------------------------------------------
% Plot bleed system with or without renewal.
figure(1)
plot(mff, ffc1, '-'); hold on
plot(mff, ffc2, '--');
plot(mff, ffmc1, 'o');
plot(mff, ffmc2, 'x'); hold off

text(41, 125, ['Usage ', num2str(usageRate), ' h/m'], 'FontSize', 13);
text(41, 75, ['POF_{ctrl} = ', num2str(pofCtrl)], 'FontSize', 13);

% Format
grid on
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');

leg = legend('w/o renewal', 'w/ renewal', 'MCS w/o renewal', 'MCS w/ renewal');
set(leg, 'Location', 'northwest');
set(leg, 'Color', 'none');
set(leg, 'FontSize', 14);
set(leg, 'FontWeight', 'bold');

xl = xlabel('Calandar Time (months)');
set(xl, 'FontSize', 15);
set(xl, 'FontWeight', 'bold');
yl = ylabel('Failure Forcast');
set(yl, 'FontSize', 15);
set(yl, 'FontWeight', 'bold');

% Plot bleed system rhythm.
figure(2)
plot(mff, [ffc1(1); ffc1(2:end) - ffc1(1:end-1)], '-'); hold on
plot(mff, [ffc2(1); ffc2(2:end) - ffc2(1:end-1)], '--');
plot(mff, [ffmc1(1); ffmc1(2:end) - ffmc1(1:end-1)], 'o');
plot(mff, [ffmc2(1); ffmc2(2:end) - ffmc2(1:end-1)], 'x'); hold off

% Format
grid on
axis([0, 60, 0, 10]);
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');

leg = legend('w/o renewal', 'w/ renewal', 'MCS w/o renewal', 'MCS w/ renewal');
set(leg, 'Location', 'southwest');
set(leg, 'Color', 'none');
set(leg, 'FontSize', 14);
set(leg, 'FontWeight', 'bold');

xl = xlabel('Calandar Time (months)');
set(xl, 'FontSize', 15);
set(xl, 'FontWeight', 'bold');
yl = ylabel('Failures Per Month');
set(yl, 'FontSize', 15);
set(yl, 'FontWeight', 'bold');