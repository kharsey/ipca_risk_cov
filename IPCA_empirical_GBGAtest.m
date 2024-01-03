
%% Set estimation parameters and data choice

% % immediately below can be commented for sbatch
% clear
% K = 1;
% bootsims = 5;
% dataname = 'IPCADATA_FNW36_RNKDMN_CON'

% Load data
load(['../Data/' dataname]);
% als_opt
als_opt.MaxIterations       = 5000;
als_opt.Tolerance           = 1e-6;


%% Start estimation

disp(['IPCA_empirical_GBGAtest starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);

[GammaBeta_XSVD,s,v]    = svds(X,K);

%% ALS
disp([' ALS started at ' datestr(clock)])

% Numerical ALS procedure: starting from QSVD guess
tic;
GB_Old      = GammaBeta_XSVD;
F_Old       = s*v';%ones(K,T);%
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
GB_GB   = GB_New;
GB_F    = F_New;
GB_L    = mean(F_New,2);
timing.als_gb.time = toc;
timing.als_gb.iter = iter;
timing.als_gb.tols = tols;
disp([' NUM_ALS for GB-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_gb.time) ' sec'])

% Numerical ALS procedure: starting from GammaBeta guess
tic;
GB_Old      = [GB_GB zeros(L,1)];
F_Old       = GB_F;
tol         = 1;
iter        = 0;
tols        = nan(500,1);
while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
    [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts,ones(1,T));
    tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
    F_Old   = F_New;
    GB_Old  = GB_New;
    iter    = iter+1;
    tols    = [tols(2:end);tol];
end
GBGA_GB      = GB_New(:,1:end-1);
GBGA_GA      = GB_New(:,end);
GBGA_F       = F_New;
GBGA_L       = mean(GBGA_F,2);
timing.als_gbga.time = toc;
timing.als_gbga.iter = iter;
timing.als_gbga.tols = tols;
disp([' NUM_ALS for GBGA-model done after ' num2str(iter) ' iterations'...
    ' at ' datestr(clock) ' after ' num2str(timing.als_gbga.time) ' sec'])


%% Fits and R2s

% Fits
RFITS_GB         = nan(N,T);
RFITS_pred_GB    = nan(N,T);
RFITS_GBGA        = nan(N,T);
RFITS_pred_GBGA   = nan(N,T);
XFITS_GBGA      = nan(L,T);
XFITS_pred_GBGA = nan(L,T);
XFITS_GB        = nan(L,T);
XFITS_pred_GB   = nan(L,T);
for t=1:T
    RFITS_GB(:,t)        = Z(:,:,t)*GB_GB*GB_F(:,t) ;
    RFITS_pred_GB(:,t)   = Z(:,:,t)*GB_GB*GB_L ;
    RFITS_GBGA(:,t)       = Z(:,:,t)*(GBGA_GA + GBGA_GB*GBGA_F(:,t) ) ;
    RFITS_pred_GBGA(:,t)  = Z(:,:,t)*(GBGA_GA + GBGA_GB*GBGA_L ) ;
    XFITS_GBGA(:,t)      = W(:,:,t)*(GBGA_GA + GBGA_GB*GBGA_F(:,t));
    XFITS_pred_GBGA(:,t) = W(:,:,t)*(GBGA_GA + GBGA_GB*GBGA_L);
    XFITS_GB(:,t)        = W(:,:,t)*(GB_GB*GB_F(:,t));
    XFITS_pred_GB(:,t)   = W(:,:,t)*(GB_GB*GB_L);
end
QFITS_GB         = GB_GB*GB_F;
QFITS_pred_GB    = repmat(GB_GB*GB_L,1,T);
QFITS_GBGA        = GBGA_GA + GBGA_GB*GBGA_F;
QFITS_pred_GBGA   = repmat(GBGA_GA + GBGA_GB*GBGA_L , 1 , T);

% R2s
xret(LOC(:)==0)         = nan;
totalsos                = mySOS(xret); 
RR2_total_GB            = 1 - mySOS( xret(LOC(:)) - RFITS_GB(LOC(:))  )/totalsos;
RR2_pred_GB             = 1 - mySOS( xret(LOC(:)) - RFITS_pred_GB(LOC(:))  )/totalsos;
RR2_total_GBGA          = 1 - mySOS( xret(LOC(:)) - RFITS_GBGA(LOC(:))  )/totalsos;
RR2_pred_GBGA           = 1 - mySOS( xret(LOC(:)) - RFITS_pred_GBGA(LOC(:))  )/totalsos;

QR2_total_GB            = 1 - mySOS( Q - QFITS_GB  )/mySOS(Q);
QR2_pred_GB             = 1 - mySOS( Q - QFITS_pred_GB  )/mySOS(Q);
QR2_total_GBGA          = 1 - mySOS( Q - QFITS_GBGA  )/mySOS(Q);
QR2_pred_GBGA           = 1 - mySOS( Q - QFITS_pred_GBGA  )/mySOS(Q);

XR2_total_GBGA       = 1 - mySOS( X - XFITS_GBGA  )/mySOS(X);
XR2_pred_GBGA        = 1 - mySOS( X - XFITS_pred_GBGA  )/mySOS(X);
XR2_total_GB         = 1 - mySOS( X - XFITS_GB  )/mySOS(X);
XR2_pred_GB          = 1 - mySOS( X - XFITS_pred_GB  )/mySOS(X);
   
disp(['  estimation completed ' datestr(clock)])


%% Bootstrap
if bootsims>0
    tic
    
    RESID = X - XFITS_GBGA;
    boot_GA = nan(bootsims,L);
    
    rng('shuffle')
    
    for boot=1:bootsims
        % bootstrap indexes
        btix    = [];
        block   = 1;        
        tmp     = 0;
        while tmp<T
            fillin = unidrnd(T) + (0:(block-1));
            fillin = fillin(fillin<=T);
            btix(tmp+(1:length(fillin))) = fillin;
            tmp = tmp+length(fillin);
        end
        btix    = btix(1:T);

        dof     = 5;
        tvar    = dof/(dof-2);

        % Data construction
        X_b     = XFITS_GB + bsxfun(@times, RESID(:,btix) , (1/sqrt(tvar))*random('t',dof,1,T) );% yes wild
        % Estimation
        GB_Old      = [GB_GB zeros(L,1)];
        F_Old       = GB_F;
        tol         = 1;
        iter        = 0;
        while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
            [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X_b,Nts,ones(1,T));
            tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
            F_Old   = F_New;
            GB_Old  = GB_New;
            iter    = iter+1;
        end
        boot_GA(boot,:) = GB_New(:,end);
    end  
    timing.boot.sims = bootsims;
    timing.boot.time = toc;
end


%% Save results
if exist('clustersuffix')==0
    clustersuffix='';
end
suffix = ['_boot' num2str(bootsims)];
save(['../Data/Results_GBGAtest_' dataname '_K' num2str(K) suffix clustersuffix] ...
    , 'xret' , 'W' , 'date' , 'LOC' , 'Nts' ...
    , 'R*' ...
    , 'Q*' ...
    , 'X*' ...
    , 'GBGA_*' , 'GB_*' ...
    , 'boot*' ...
    , 'timing' ); 





































