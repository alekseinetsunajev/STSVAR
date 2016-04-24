clear
% uncomment next line to use parallel computing
%matlabpool open;
Nrep = 20;                  % number of bootstrap replications
q = 0.84;                     % 68% confidence bands

%% VAR part
% obtain matrices of restrictions
[R_Matr_SR_BLeit, R_Matr_LR, R_Matr_SR_RC] = get_r_matr;

% estimate SVAR using Bjornland-Leitemo restrictions
VAR_BLeitemo =  classVAR(1, 1,  R_Matr_SR_BLeit,  R_Matr_LR);
VAR_BLeitemo = VAR_BLeitemo.estimateSVAR;

% bootstrap loop
for ii = 1: Nrep
    VAR_BLeitemo_BIRS(ii) = VAR_BLeitemo;
end

parfor ii = 1: Nrep
    disp(ii)
    VAR_BLeitemo_BIRS(ii) = VAR_BLeitemo_BIRS(ii).mbb_bootstrap;
end

for ii = 1:Nrep
    BIRS_SR_BL(:,:, ii) = VAR_BLeitemo_BIRS(ii).IR_SRun - VAR_BLeitemo.IR_SRun;
    BIRS_LR_BL(:,:, ii) = VAR_BLeitemo_BIRS(ii).IR_LRun - VAR_BLeitemo.IR_LRun;
end

% display results
quantH = (quantile(BIRS_SR_BL,q, 3));
quantL = (quantile(BIRS_SR_BL,1-q, 3));
quantH_LR = (quantile(BIRS_LR_BL,q, 3));
quantL_LR = (quantile(BIRS_LR_BL,1-q, 3));

plot_IR = VAR_BLeitemo.IR_SRun;
plot_IR([4:5:25], :) = VAR_BLeitemo.IR_LRun([4:5:25], :);

Q_L = VAR_BLeitemo.IR_SRun - quantH;
Q_L([4:5:25], :) = VAR_BLeitemo.IR_LRun([4:5:25], :) - quantH_LR([4:5:25], :);

Q_H = VAR_BLeitemo.IR_SRun - quantL;
Q_H([4:5:25], :) = VAR_BLeitemo.IR_LRun([4:5:25], :) - quantL_LR([4:5:25], :);

plot_ir_v(plot_IR, Q_L, Q_H, [1 2 3 4 5], VAR_BLeitemo.h, VAR_BLeitemo.T, 1)
