% Data pre-processing for case study 2 of Weibull Handbook

clear; clc;

% Bases     A  B  C  A  A  C  C  E  F
baseCode = [1, 2, 3, 1, 1, 3, 3, 5, 6];
baseCode = [baseCode, 4*ones(1, 10)];

% Time of failure
tf = [153,872, 1568, 1428, 212, 808, 1405, 64, 32, ...
    1198, 884, 1251, 1249, 708, 1082, 884, 1105, 828, 1013];

% Non-failure time
tt = 250:100:2150;
nnf = [2, 0, 0, 2, 2, 9, 23, 27, 20, 22, 22, 11, 11, 20, 8, 4, 2, 3, 3, 1];

% Assemble data for Weibull analysis. -------------------------------------
% Indices for data without base D
idxDataWoD = baseCode ~= 4;

% Time of failure without base D
tfWoD = tf(idxDataWoD)';
idxFailedWoD = ones(size(tfWoD));

% Time of failure for base D only
tfForD = tf(~idxDataWoD)';
idxFailedForD = ones(size(tfForD));

% Time of non-failure at base D
tnfForD = nan(sum(nnf), 1);
for i = 1:numel(tt)
    if nnf(i) ~= 0
        idxStart = sum(~isnan(tnfForD)) + 1;
        idxEnd = idxStart + nnf(i) - 1;
        
        tnfForD(idxStart:idxEnd) = tt(i);
    end
end

% Data without base D
tWoD = tfWoD;
%tWoD = [tfWoD; tnfForD];
%idxFailedWoD = [idxFailedWoD; zeros(size(tnfForD))];

% Data for base D
tForD = [tfForD; tnfForD];
idxFailedForD = [idxFailedForD; zeros(size(tnfForD))];

% Save data for Weibull analysis.
csvwrite('case_study_2_data_wo_d.csv', [tWoD, idxFailedWoD]);
csvwrite('case_study_2_data_for_d.csv', [tForD, idxFailedForD]);

% Plot for data visualization.
bar(tt, nnf);
histogram([tfWoD; tForD], 'Normalization', 'probability');

% Format
grid on
set(gca, 'FontSize', 14);
set(gca, 'FontWeight', 'bold');

xl = xlabel('TSN (h)');
set(xl, 'FontSize', 15);
set(xl, 'FontWeight', 'bold');
yl = ylabel('Probability');
set(yl, 'FontSize', 15);
set(yl, 'FontWeight', 'bold');