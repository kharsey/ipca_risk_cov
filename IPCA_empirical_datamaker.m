% IPCA_empirical_datamaker
clear


%% mat file name
% the specification choices in this section with automatically add suffixes to this name 
mat_file_name = 'IPCADATA_FNW36'

datacallname = 'IPCA_empirical_datacall_FNW36'

option_rnkdmn_beg                   = 1; % find unit ranks before anything else [RNKDMN in name]
option_zscore_beg                   = 0; % makes every char into a CS zscore [ZSCORE in name]

% only 1 of following 2 can be turned on
option_tsmean                       = 0; % calculate the tsmean [TS in name]
option_histmean                     = 0; % calculate the histmean [HIST in name]

% only 1 of following 3 can be turned on
option_mean_dataset                 = 0; % just the means [MEAN in name]
option_dev_dataset                  = 0; % just the devs [DEV in name]
option_meandev_dataset              = 0; % means and devs [MEANDEV in name]

option_constant                     = 1; % include a constant characteristic at the end
option_orthogonalize                = 0; % orthogonalize all chars every period

option_keepthese                    = {};
% {'beta' 'lme' 'at' 'beme' 'e2p' 'prof' 'idio_vol' 'roe' 'investment' 'cum_return_12_2' 'cum_return_36_13' 'cum_return_1_0' 'rel_to_high_price'};
%{'beta' 'lme' 'at' 'cum_return_1_0' 'cum_return_12_2' 'lturnover' ...
% ...%                                         'rel_to_high_price' 'cum_return_36_13' 'suv' 'idio_vol'};
%     {...
%     'beta'    'e2p'    'beme'    'q'    'a2me'    'lme'    'cto'    'oa'    'roe' 'lturnover'    'noa'    ...
%     'idio_vol'    'ato'    'dpi2a'    'investment'    'pm'    'rna'    'suv'    'roa'    'free_cf'    'ol'    'c' 'spread_mean'    'at'    'lev'    'prof'    's2p' ...
%     'sga2m'    'd2a'    'fc2y'    'pcm' 'rel_to_high_price' 'cum_return_36_13' };%'cum_return_1_0'};% 'cum_return_12_2'};
%     {'beta' 'cum_return_1_0' 'lme' 'cum_return_12_2' 'rel_to_high_price' 'cum_return_36_13' 'cum_return_12_7' 'at'};
option_keepbig                      = 0;
option_keepsmall                    = 0;
option_keepbighalf                  = 0;
option_keepsmallhalf                = 0;
option_subsample                    = [];% either empty or 2-vector with yyyymm specified range (eg [198001 199912])
option_keepfirsthalf                = 0;
option_keepsecondhalf               = 0;
option_keeprandomhalf               = 0;% 0 means no; 1 means use random sample vector; 2 means complement of random sample vector

%% Datacall
% Whatever script is called needs to produce two arrays: ret and chars
% xret needs to be NxT
% chars needs to be NxTxL
% chars needs to be lagged: if col j of ret is the monthly excess return at the end of month t, col j of chars is the
% characteristics as of the end of month t-1

disp(['IPCA_empirical_datamaker started ' datestr(clock)])

run(datacallname);

[N,T,L] = size(chars);
disp(['  data read in ' datestr(clock)])

% find stocks each month that have all the chars and rets
LOC = ~isnan(chars);
LOC = all(LOC,3)&~isnan(xret);
disp(['  LOC found ' datestr(clock)])

% Check for a big enough cross section
keepthese = sum(LOC)>=max(100,(option_meandev_dataset+1)*L+option_constant+1);
LOC     = LOC(:,keepthese);
chars   = chars(:,keepthese,:);
xret    = xret(:,keepthese);
date    = date(keepthese);
[N,T,L] = size(chars);



%% keepthese
if ~isempty(option_keepthese)
    loc = [];
    for j=1:length(option_keepthese)
        loc = [loc find(strcmp(option_keepthese{j},charnames))];
    end
    chars       = chars(:,:,loc);
    [N,T,L]     = size(chars);
    charnames = charnames(loc);
    disp(['We decided to only keep ' charnames])
end


%% transform chars to RNKDMN at the beginning
if option_rnkdmn_beg
chars = tiedrank(chars);
chars = (chars-1)./(max(chars)-1);
chars = chars-0.5;
disp(['  You found RNKDMN at the beginning at ' datestr(clock)])
mat_file_name = [mat_file_name '_RNKDMN'];
end


%% transform chars to ZSCORE at the beginning
if option_zscore_beg
chars = zscore2(chars);
disp(['  You found ZSCORE at the beginning at ' datestr(clock)])
mat_file_name = [mat_file_name '_ZSCORE'];
end

%% FIND TSMEAN or HISTMEAN
if option_tsmean && option_histmean
    error('option_tsmean and option_histmean are both turned on')
end

if option_tsmean
tmp = repmat(squeeze(nanmean(chars,2)),1,1,T);% tmp is NxLxT
tmp = permute(tmp,[1 3 2]);% tmp is NxTxL
tmp(isnan(chars(:))) = nan;
disp('  You found TS')
meansuffix = '_TS';
end

if option_histmean
tmp               = zeros(size(chars));
tmploc            = repmat(LOC,1,1,L);
tmp(tmploc(:))    = chars(tmploc(:));
tmp2              = cumsum(tmp)./cumsum(tmploc);
tmp2(~tmploc(:))  = nan;
tmp               = tmp2;
disp('  You found HIST')
meansuffix = '_HIST';
end

%% MAKE a MEAN or DEV or MEANDEV or DEVMEAN data set
if option_mean_dataset + option_dev_dataset + option_meandev_dataset > 1
    error('more than one of mean/dev/meandev are turned on')
end

if option_mean_dataset
chars = tmp;
clear tmp*
disp('  You found MEAN')
mat_file_name = [mat_file_name meansuffix '_MEAN'];
end

if option_dev_dataset
tmp2                = chars-tmp;
chars               = tmp2;
clear tmp*
disp('  You found DEV')
mat_file_name = [mat_file_name meansuffix '_DEV'];
end

if option_meandev_dataset
tmp2                = nan(N,T,2*L);
tmp2(:,:,1:L)       = tmp;
tmp2(:,:,L+1:end)   = chars-tmp;
chars               = tmp2;
L = 2*L;
clear tmp*
disp('  You found MEANDEV')
mat_file_name = [mat_file_name meansuffix '_MEANDEV'];
end





%% add a constant as the last characteristic
if option_constant
[N,T,L]         = size(chars);
chars(:,:,L+1)  = 1;
[N,T,L]         = size(chars);
disp(['  You added a constant at ' datestr(clock)])
mat_file_name = [mat_file_name '_CON'];
end


%% KEEPBIG or KEEPSMALL
if option_keepbig || option_keepsmall
    chars = permute(chars,[1 3 2]);% make chars NxLxT
    if option_keepbig && option_keepsmall
        disp(['Problem: both keepbig and keepsmall were selected -- going with "keepbig"'])
        option_keepsmall = 0;
    end
    small1000 = LOC;
    for t=1:T
        tmp             = sort(chars(LOC(:,t),8,t)); %% REQUIRES SIZE TO BE THE 8th char!!!
        tmp2            = tmp(max(1,length(tmp)-999));
        small1000(:,t)  = small1000(:,t) & chars(:,8,t)<=tmp2;
    end
    big1000 = LOC & ~small1000;
    chars = permute(chars,[1 3 2]);% make chars NxTxL
    if option_keepbig
        LOC = big1000;
        disp(['  You keep the big stocks at ' datestr(clock)])
        mat_file_name = [mat_file_name '_KEEPBIG'];
    else
        LOC = small1000;
        % Check for a big enough cross section
        keepthese = sum(LOC)>=max(100,(option_meandev_dataset+1)*L+option_constant+1);
        LOC     = LOC(:,keepthese);
        chars   = chars(:,keepthese,:);
        xret    = xret(:,keepthese);
        date    = date(keepthese);
        [N,T,L] = size(chars);
        disp(['  You keep the small stocks at ' datestr(clock)])
        mat_file_name = [mat_file_name '_KEEPSMALL'];
    end
end


%% KEEPBIGHALF or KEEPSMALLHALF
if option_keepbighalf || option_keepsmallhalf
    chars = permute(chars,[1 3 2]);% make chars NxLxT
    if option_keepbighalf && option_keepsmallhalf
        disp(['Problem: both keepbighalf and keepsmallhalf were selected -- going with "keepbighalf"'])
        option_keepsmallhalf = 0;
    end
    smallhalf = LOC;
    for t=1:T
        tmp             = sort(chars(LOC(:,t),8,t)); %% REQUIRES SIZE TO BE THE 8th char!!!
        tmp2            = tmp(max(1,length(tmp)-999));
        smallhalf(:,t)  = smallhalf(:,t) & chars(:,8,t) < nanmedian(chars(:,8,t));
    end
    bighalf = LOC & ~smallhalf;
    chars = permute(chars,[1 3 2]);% make chars NxTxL
    if option_keepbighalf
        LOC = bighalf;
        disp(['  You keep the big stocks at ' datestr(clock)])
        mat_file_name = [mat_file_name '_KEEPBIGHALF'];
    else
        LOC = smallhalf;
        % Check for a big enough cross section
        keepthese = sum(LOC)>=max(100,(option_meandev_dataset+1)*L+option_constant+1);
        LOC     = LOC(:,keepthese);
        chars   = chars(:,keepthese,:);
        xret    = xret(:,keepthese);
        date    = date(keepthese);
        [N,T,L] = size(chars);
        disp(['  You keep the small stocks at ' datestr(clock)])
        mat_file_name = [mat_file_name '_KEEPSMALLHALF'];
    end
end

%% SUBSAMPLE
if ~isempty(option_subsample)
    t1 = find(date>=option_subsample(1),1,'first');
    t2 = find(date<=option_subsample(2),1,'last');
    option_subsample(1) = date(t1);
    option_subsample(2) = date(t2);
    LOC = LOC(:,t1:t2);
    chars = chars(:,t1:t2,:);
    [N,T,L] = size(chars);
    xret = xret(:,t1:t2);
    date = date(t1:t2);
    mat_file_name = [mat_file_name '_' num2str(option_subsample(1)) '-' num2str(option_subsample(2))];
end

%% KEEPFIRSTHALF or KEEPSECONDHALF
if option_keepfirsthalf || option_keepsecondhalf
    loc     = find(LOC(:));
    loc     = loc(floor(length(loc)/2));
    loc     = ceil(loc/N);
    if option_keepfirsthalf
        LOC = LOC(:,1:loc);
        chars = chars(:,1:loc,:);
        xret = xret(:,1:loc);
        date = date(1:loc);
        [N,T,L] = size(chars);
        mat_file_name = [mat_file_name '_KEEPFIRSTHALF'];
    elseif option_keepsecondhalf
        LOC = LOC(:,loc+1:end);
        chars = chars(:,loc+1:end,:);
        xret = xret(:,loc+1:end);
        date = date(loc+1:end);
        [N,T,L] = size(chars);
        mat_file_name = [mat_file_name '_KEEPSECONDHALF'];
    end
end


%% KEEPRANDOMHALF
if option_keeprandomhalf
    rh = load(['../Data/randomhalf_' mat_file_name]);
    switch option_keeprandomhalf
        case 1
            LOC(rh.randomhalfcomp,:) = false;
        case 2
            LOC(rh.randomhalf,:) = false;
        otherwise
            error('What option was passed?')
    end
    clear rh
    mat_file_name = [mat_file_name '_RANDOMHALF' num2str(option_keeprandomhalf)];
end
    



%% Orthogonalize or not

Nts = sum(LOC);

% ORTHOGONALIZATION
if option_orthogonalize
    % orthogonalize chars
    disp(['  starting to orthogonalize Z ' datestr(clock)])
    Z = nan(T,N,L);% Z is TxNxL
    for t=1:T %parfor makes it too big
        slice = squeeze( chars(LOC(:,t),t,:) );
        if sum(isnan(slice(:)))>0; error('wtf');end

        % Cholesky orthogonalize
        slice = slice/chol(slice'*slice);
        slice = slice*diag(sign(mean(slice)));

        slice2 = nan(N,L);
        slice2(LOC(:,t),:) = slice; 
        Z(t,:,:) = slice2;
        Z(t,:,:) = Z(t,:,:) * sqrt(Nts(t)); % we want (1/N)Zt'*Zt = I_L
    end
    Z       = permute(Z,[2 3 1]); % Z is now NxLxT
    chars   = permute(chars,[1 3 2]); % chars is now NxLxT
    disp(['  Z orthogonalized ' datestr(clock)])

    % Construct X
    Q = nan(L,T);
    for t=1:T % parfor makes it too big
       Q(:,t) = (1/Nts(t)) * Z(LOC(:,t),:,t)'*xret(LOC(:,t),t) ;
    end
    X       = Q;
    W       = repmat(eye(L),1,1,T);
    mat_file_name = [mat_file_name '_ORTHO'];
% NO ORTHOGONALIZATON    
else
    chars   = permute(chars,[1 3 2]); % chars is now NxLxT
    Z       = chars;
    % Construct X, W, Xtil
    W       = nan(L,L,T);
    Q       = nan(L,T);
    X       = nan(L,T);
    for t=1:T % parfor makes it too big
        W(:,:,t)    = (1/Nts(t)) * Z(LOC(:,t),:,t)' * Z(LOC(:,t),:,t);
        X(:,t)   = (1/Nts(t)) * Z(LOC(:,t),:,t)' * xret(LOC(:,t),t);
        Q(:,t)      = W(:,:,t) \ X(:,t) ;
    end
end
      

%% Save mat file
disp(['Saving ' mat_file_name])
save(['../Data/' mat_file_name],'-v7.3')
disp(['IPCA_empirical datamaker done at ' datestr(clock)])

