% Figure 4-7. Bleed System Risk - Base D

clear; clc;

% Constant
usageRate = 25;    % Operation time per month, [h/m]

% Load data.
data = csvread('case_study_2_data_for_d.csv');
tsn = data(:, 1);
idxFailed = data(:, 2);

% Replace failed part w/ new one.
%tsn(logical(idxFailed)) = 0;

% Consider only suspension data.
tsn(logical(idxFailed)) = [];

% Parameters
etaParam = 2004.45;
betaParam = 5.239307;

% Model for failure forcast
FF = @(t0, t1) (wblcdf(t1, etaParam, betaParam) - ...
    wblcdf(t0, etaParam, betaParam))./(1 - wblcdf(t0, etaParam, betaParam));

% Months for failure forcast
mff = 1:60;
tff = mff*usageRate;

% Do calculation.
ff = zeros(numel(mff), 1);
for i = 1:numel(mff)
    ff(i) = sum(FF(tsn, tsn + tff(i)));
end

% Plot Figure 4-7.
plot(mff, ff, '-^');

text(41, 75, ['Usage ', num2str(usageRate), ' h/m'], 'FontSize', 13);

% Format
grid on
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');

xl = xlabel('Calandar Time (months)');
set(xl, 'FontSize', 15);
set(xl, 'FontWeight', 'bold');
yl = ylabel('Failure Forcast');
set(yl, 'FontSize', 15);
set(yl, 'FontWeight', 'bold');