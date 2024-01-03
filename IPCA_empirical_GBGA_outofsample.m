
%% Set estimation parameters and data choice

% % the following can be commented out if using sbatch
% clear
% dataname = 'IPCADATA_FNW36_RNKDMN_CON'
% for K=1:6
%     clearvars -except K dataname


% als_opt
als_opt.MaxIterations       = 5000;
als_opt.Tolerance           = 1e-6;
startindex = 60;


%% Start estimation

disp(['IPCA_empirical_GB starting at ' datestr(clock) ': K=' num2str(K) ', data=' dataname]);

% Load data
load(['../Data/' dataname]);


bigX    = X;
bigW    = W;
bigNts  = Nts;

%% Estimation

OOSRFITS_pred_GB    = nan(N,T);
OOSXFITS_pred_GB = nan(L,T);
OOSRealFact         = nan(K,T);
OOSRealTan          = nan(T,1);
OOSRFITS_GB         = nan(N,T);
OOSXFITS_GB      = nan(L,T);
OOSARBPTF           = nan(1,T);


for t=startindex:T-1
    X       = bigX(:,1:t);% this is X known through t
    W       = bigW(:,:,1:t);% this is W known through t
    Nts     = bigNts(1:t);% this is Nts known through t

    [GammaBeta_XSVD,s,v]    = svds(X,K);
    % Numerical ALS procedure -- GB
    tic;
    if t==startindex
        GB_Old      = GammaBeta_XSVD;
        F_Old       = s*v';%ones(K,T);%
    else
        GB_Old      = GammaBeta;
        F_Old       = [FactorGB s*v(end,:)'];
    end
    tol         = 1;
    iter        = 0;
    while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
        [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts);
        tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
        F_Old   = F_New;
        GB_Old  = GB_New;
        iter    = iter+1;
    end
    FactorGB    = F_New;
    GammaBetaGB = GB_New;
    LambdaGB    = mean(F_New,2);
    disp([' NUM_ALS for GB-model for t=' num2str(t) ' done after ' num2str(iter) ' iterations'...
        ' at ' datestr(clock) ' after ' num2str(toc) ' sec'])
    % Numerical ALS procedure -- GBGA
    tic;
    if t==startindex
        GB_Old      = [GammaBeta_XSVD zeros(L,1)];
        F_Old       = s*v';%ones(K,T);%
    else
        GB_Old      = [GammaBeta GammaAlpha];
        F_Old       = [Factor s*v(end,:)'];
    end
    tol         = 1;
    iter        = 0;
    while iter<=als_opt.MaxIterations && tol>als_opt.Tolerance
        [GB_New,F_New] = num_IPCA_estimate_ALS(GB_Old,W,X,Nts,ones(1,t));
        tol     = max([ abs(GB_New(:)-GB_Old(:)) ; abs(F_New(:)-F_Old(:)) ]);
        F_Old   = F_New;
        GB_Old  = GB_New;
        iter    = iter+1;
    end
    GammaBeta  = GB_New(:,1:end-1);
    GammaAlpha = GB_New(:,end);
    Factor     = F_New;
    Lambda     = mean(F_New,2);
    disp([' NUM_ALS for GBGA-model for t=' num2str(t) ' done after ' num2str(iter) ' iterations'...
        ' at ' datestr(clock) ' after ' num2str(toc) ' sec'])
    
    % Arb ptf
    tmp = ( GammaAlpha'/(Z(LOC(:,t+1),:,t+1)'*Z(LOC(:,t+1),:,t+1))*Z(LOC(:,t+1),:,t+1)' );
    tmp2 = xret(LOC(:,t+1),t+1);
    OOSARBPTF(t+1)    = tmp*tmp2;

    % Tan using GB factors and arb
    % training sample arb ptf returns
    ts_arbptf = nan(1,t);
    for tt = 1:t
        tmp = ( GammaAlpha'/(Z(LOC(:,tt),:,tt)'*Z(LOC(:,tt),:,tt))*Z(LOC(:,tt),:,tt)' );
        tmp2 = xret(LOC(:,tt),tt);
        ts_arbptf(tt) = tmp*tmp2;
    end
    OOSReal_FGB(:,t+1) = ( GammaBetaGB'*bigW(:,:,t+1)*GammaBetaGB )\( GammaBetaGB'*bigX(:,t+1) );
    OOSReal_FGB_Arb_Tan(t+1) = tanptfnext([FactorGB; ts_arbptf]', [OOSReal_FGB(:,t+1); OOSARBPTF(t+1)]');
    OOSReal_FGB_Tan(t+1) = tanptfnext(FactorGB', OOSReal_FGB(:,t+1)');
end


X = bigX;
W = bigW;
Nts = bigNts;

disp(['  estimation completed ' datestr(clock)])


%% Save results
save(['../Data/Results_GBGA_outofsample_' dataname '_K' num2str(K)] ...
    , 'xret' , 'W' , 'X' , 'date' , 'LOC' , 'Nts' ...
    , 'OOS*' ); 


% end


































