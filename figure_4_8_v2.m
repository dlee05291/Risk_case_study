% Figure 4-8. Bleed System With Renewal
%
% Revision history
% 022020 LDY Code was created. Renewal was accomplished by adding new units
%            to suspension fleet when POF of individual unit exceed control
%            level.
% 022120 LDY Code was modified to match results from conventional failure
%            forcast method with results from MCS.

clear; clc;

% Constant
usageRate = 50;    % Operation time per month, [h/m]
nsample = 1e6;     % Number of samples for MCS

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
nfailed = zeros(numel(tff), 1);

% Do calculation.
for i = 1:numel(tff)
    nfailed(i) = sum(FF(tsn, tsn + tff(i)));
end

% Get failure forecast.
ffc1 = nfailed;

% With renewal (conventional) ---------------------------------------------
% Initialize variables for the number of failed unit.
nfailed1 = zeros(numel(tff), 1);
deltaffmc2 = zeros(numel(tff), 1);
deltaNew = zeros(numel(tff), 1);

% Do calculation.
for i = 1:numel(tff)
    % Number of failure for initial unit
    nfailed1(i) = sum(FF(tsn, tsn + tff(i)));
    
    % Number of failure for renewed unit
    if i == 1
        % Number of new units introduced
        deltaNew(i) = nfailed1(i);
    else
        % Number of new units introduced
        deltaNew(i) = nfailed1(i) - nfailed1(i - 1);
        
        % Increase in probability of failure for current step
        deltaPOF = FF(0, flip(tff(1:i-1)));
        
        % Increase in number of failure for current step.
        deltaffmc2(i) = sum(deltaPOF(:).*deltaNew(1:i-1));
    end
end

% Get failure forecast.
ffc2 = nfailed1 + deltaffmc2;

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
pof2 = zeros(numel(tff), 1);

% Do MCS.
% Time to failure
ttf2 = TTF(nsample);

for i = 1:numel(tff)
    % Indices for failed unit
    idxFailedPerUnit = ttf2 < tff(i);
    
    % Get POF for renewed unit.
    pof2(i) = sum(idxFailedPerUnit)/nsample;
    
    % Number of failure for renewed unit
    if i == 1
        % Number of new units introduced
        deltaNew(i) = ffmc1(i);
    else
        % Number of new units introduced
        deltaNew(i) = ffmc1(i) - ffmc1(i - 1);
        
        % Increase in probability of failure for current step
        deltaPOF = flip(pof2(1:i-1));
        
        % Increase in number of failure for current step.
        deltaffmc2(i) = sum(deltaPOF(:).*deltaNew(1:i-1));
    end
end

% Get failure forecast.
ffmc2 = ffmc1 + deltaffmc2;

% Plot --------------------------------------------------------------------
% Plot bleed system with or without renewal.
figure(1)
plot(mff, ffc1, '-'); hold on
plot(mff, ffc2, '--');
plot(mff, ffmc1, 'o');
plot(mff, ffmc2, 'x'); hold off

text(41, 125, ['Usage ', num2str(usageRate), ' h/m'], 'FontSize', 13);

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

% Plot number of units.
mff2 = [0, mff];

figure(3)
h(1) = plot(mff2, nunit - [0; nfailed1] + cumsum([0; deltaNew]), '-'); hold on
h(2) = plot(mff2, nunit - [0; nfailed1], '--');
h(3) = plot(mff2, cumsum([0; deltaNew]), '-.'); hold off

% Format
grid on
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
set(h, 'LineWidth', 2);

leg = legend('w/ renewal', 'w/o renewal', 'Number of renewal');
set(leg, 'Location', 'best');
set(leg, 'Color', 'none');
set(leg, 'FontSize', 14);
set(leg, 'FontWeight', 'bold');

xl = xlabel('Calandar Time (months)');
set(xl, 'FontSize', 15);
set(xl, 'FontWeight', 'bold');
yl = ylabel('Number of Units');
set(yl, 'FontSize', 15);
set(yl, 'FontWeight', 'bold');