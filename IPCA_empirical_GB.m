
%% Set estimation parameters and data choice

Krange = [1:6];

for K=Krange
    K
    
if K==Krange(1)
    dataname = 'IPCADATA_FNW36_RNKDMN_CON'
    % Load data
    clearvars -except Krange K dataname
    load(['../Data/' dataname]);
    % als_opt
    als_opt.MaxIterations       = 10000;
    als_opt.Tolerance           = 1e-6;
end

%% Start estimation

disp(['IPCA_empirical_GB starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);

[GammaBeta_initial,s,v]    = svds(X,K);

%% ALS
disp([' ALS started at ' datestr(clock)])

% Numerical ALS procedure: starting from QSVD guess or previous factor guess
tic;
if exist('Factor')==1
    GB_Old  = [GammaBeta GammaBeta_initial(:,K)];
    F_Old   = [Factor; ones(1,T)];
else
    GB_Old      = GammaBeta_initial;
    F_Old       = s*v';%ones(K,T);%
end
tol         = 1;
iter        = 0;
tols        = nan(500,1);
while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
    [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts);
    tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
    F_Old   = F_New;
    GB_Old  = GB_New;
    iter    = iter+1;
    tols    = [tols(2:end);tol];
end
GammaBeta  = GB_New;
Factor     = F_New;
Lambda     = mean(F_New,2);
timing.als_xsvd.time = toc;
timing.als_xsvd.iter = iter;
timing.als_xsvd.tols = tols;
disp([' NUM_ALS for GB-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_xsvd.time) ' sec'])


%% Fits and R2s

% Fits
RFITS_GB            = nan(N,T);
RFITS_pred_GB       = nan(N,T);
XFITS_GB            = nan(L,T);
XFITS_pred_GB       = nan(L,T);
for t=1:T
    RFITS_GB(:,t)        = Z(:,:,t)*GammaBeta*Factor(:,t);
    RFITS_pred_GB(:,t)   = Z(:,:,t)*GammaBeta*Lambda;
    XFITS_GB(:,t)     = W(:,:,t)*GammaBeta*Factor(:,t);
    XFITS_pred_GB(:,t)   = W(:,:,t)*GammaBeta*Lambda;
end

% R2s
xret(LOC(:)==0)         = nan;
totalsos                = mySOS(xret); 
RR2_total_GB             = 1 - mySOS( xret(LOC(:)) - RFITS_GB(LOC(:))  )/totalsos;
RR2_pred_GB              = 1 - mySOS( xret(LOC(:)) - RFITS_pred_GB(LOC(:))  )/totalsos;

XR2_total_GB            = 1 - mySOS( X - XFITS_GB  )/mySOS(X);
XR2_pred_GB             = 1 - mySOS( X - XFITS_pred_GB  )/mySOS(X);
   
disp(['  estimation completed ' datestr(clock)])


%% Save results
save(['../Data/Results_GB_' dataname '_K' num2str(K)] ...
    , 'xret' , 'W' , 'date' , 'LOC' , 'Nts' ...
    , 'Gamma*' , 'Factor*' , 'Lambda*' ...
    , 'R*' ...
    , 'X*' ...
    , 'timing' ); 

disp('XR2_total_GB XR2_pred_GB')
disp(num2str([XR2_total_GB XR2_pred_GB]))
disp('RR2_total_GB RR2_pred_GB')
disp(num2str([RR2_total_GB RR2_pred_GB]))


end



































