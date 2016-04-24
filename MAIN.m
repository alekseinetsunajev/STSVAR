%% This is the code for the paper Structural Vector Autoregressions with Smooth Transition in Variances 
% The Interaction Between U.S. Monetary Policy and the Stock Market
% by Helmut Luetkepohl and Aleksei Netðunajev
% Code by Aleksei Netsunajev
%% general part
clear
% Uncomment next line to use parallel computing for bootstrap
% matlabpool open;
Nrep = 1000;                  % number of bootstrap replications
q = 0.84;                     % 68% confidence bands for impulse responses

%% STAR part
% Restrictions (if needed)
[R_Matr_SR_BLeit, R_Matr_LR, R_Matr_SR_RC] = get_r_matr;

TR_var_time = (1:1:447)';   % time as transition variable
TR_var_inf = xlsread('TR_var_Inf_standardized2lags.xls');       % inflation as transition variable

l_i = 0;
l_j = 0;

% Initialize objects to optimize
% These are the STSVAR models without restrictions, i.e. identified only via
% changes in volatility.
STSVAR_T = classSTSVAR(-2.77, 167, 0, 0, 0, [], [], l_i, l_j, TR_var_time);
STSVAR_INF = classSTSVAR(0.49, 4.41, 0, 0, 0, [], [], l_i, l_j, TR_var_inf);

% perform the optimization of the ST SVAR models
STSVAR_T = STSVAR_T.optimize;
STSVAR_INF = STSVAR_INF.optimize;

% compute std.errors and IRs 
STSVAR_T = STSVAR_T.std_err;
STSVAR_T = STSVAR_T.IR;

STSVAR_INF = STSVAR_INF.std_err;
STSVAR_INF = STSVAR_INF.IR;

% plot some useful information
STSVAR_T.plot_tr;
STSVAR_INF.plot_tr;
STSVAR_T.plot_u;
STSVAR_INF.plot_u;

%% bootstrap loop 
% SELECT WHAT MODEL YOU WANT TO COMPUTE IMPULSE RESPONSES YOURSELF
STSVAR = STSVAR_T;                      % comment/uncomment to use further model with TIME as transition variable
%STSVAR = STSVAR_INF;               % comment/uncomment to use further model with INFLATION(lag 2) as transition variable

% intialize STSVAR objects for bootstrapping
for jj = 1: Nrep
    BIRS(jj) = STSVAR;
end

% Bootstrapping
parfor jj = 1: Nrep
    disp(jj)
    BIRS(jj) = BIRS(jj).bootstr;
end

% plot IRS
for ii = 1 :Nrep
    B_IR_SR(:,:, ii) = BIRS(ii).IR_SRun;
    B_IR_LR(:,:, ii) = BIRS(ii).IR_LRun;
end

for ii = 1 :Nrep
    B_IR_SR(:,:, ii) = B_IR_SR(:,:, ii) - STSVAR.IR_SRun;
    B_IR_LR(:,:, ii) = B_IR_LR(:,:, ii) - STSVAR.IR_LRun;
end

q = 0.84;

quantH = (quantile(B_IR_SR,q, 3));
quantL = (quantile(B_IR_SR,1-q, 3));
quantH_LR = (quantile(B_IR_LR,q, 3));
quantL_LR = (quantile(B_IR_LR,1-q, 3));

plot_IR = STSVAR.IR_SRun;
plot_IR([4:5:25], :) = STSVAR.IR_LRun([4:5:25], :);

% the responses for stock returns should be accumulated
Q_L = STSVAR.IR_SRun - quantH;
Q_L([4:5:25], :) = STSVAR.IR_LRun([4:5:25], :) - quantH_LR([4:5:25], :);

Q_H = STSVAR.IR_SRun - quantL;
Q_H([4:5:25], :) = STSVAR.IR_LRun([4:5:25], :) - quantL_LR([4:5:25], :);

% plot function
plot_ir_v(plot_IR, Q_L, Q_H, [1 2 3 4 5], STSVAR.h, STSVAR.T, 1)
