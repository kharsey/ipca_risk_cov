
%% Set estimation parameters and data choice
clear

% the particular dataname only matters for 'QX' mode -- otherwise, it is used simply to
% bring in the xret, LOC, and date variables; for 'QX' mode, the K choice matters 
dataname = 'Results_GB_IPCADATA_FNW36_RNKDMN_CON_K1'; 
oosflag = false;
QXorR = 'R';% 'Q', 'R', 'X'
annualize = false; % annualize FF -- assumes "date" vec specifies first of the twelve months over which return calc'd
ObsChoice = 'FF';% 'FF', 'SY', 'HXZ', 'BS'

switch ObsChoice
    case 'FF'
    %% Load data
    load(['../Data/' dataname]);

    if annualize
        ffsuffix = '_ANNRET';
    else
        ffsuffix = '';
    end
    ffdata = load(['../Data/F-F_Research_Data_5_Factors_2x3_plusMOM' ffsuffix]);



    %% Start estimation
    disp(['Observable_Factor_Regressions starting at ' datestr(clock) ': data=' dataname]);



    %% FF models
    [date,loc1,loc2] = intersect(ffdata.dates,date);


    FF1 = [ffdata.Mkt_RF(loc1)];
    FF3 = [ffdata.Mkt_RF(loc1) ffdata.SMB(loc1) ffdata.HML(loc1)];
    FF4 = [FF3 ffdata.MOM(loc1)];
    FF5 = [FF3 ffdata.RMW(loc1) ffdata.CMA(loc1)];
    FF6 = [FF5 ffdata.MOM(loc1)];


    switch QXorR
        case 'R'
            RET     = xret(:,loc2);
            Rlogic  = true;
        case 'Q'
            RET     = Q(:,loc2);
            Rlogic  = false;
        case 'X'
            RET     = X(:,loc2);
            Rlogic  = false;
    end

    % preallocate
    FITS_FF1 = nan(size(RET));
    FITS_FF3 = nan(size(RET));
    FITS_FF4 = nan(size(RET));
    FITS_FF5 = nan(size(RET));
    FITS_FF6 = nan(size(RET));

    FITS_cond_FF1 = nan(size(RET));
    FITS_cond_FF3 = nan(size(RET));
    FITS_cond_FF4 = nan(size(RET));
    FITS_cond_FF5 = nan(size(RET));
    FITS_cond_FF6 = nan(size(RET));

    FF1_mean = repmat(mean(FF1),size(RET,2),1);
    FF3_mean = repmat(mean(FF3),size(RET,2),1);
    FF4_mean = repmat(mean(FF4),size(RET,2),1);
    FF5_mean = repmat(mean(FF5),size(RET,2),1);
    FF6_mean = repmat(mean(FF6),size(RET,2),1);

    FITS_OOSTan_FF1 = nan(size(RET,2),1);
    FITS_OOSTan_FF2 = nan(size(RET,2),1);
    FITS_OOSTan_FF3 = nan(size(RET,2),1);
    FITS_OOSTan_FF4 = nan(size(RET,2),1);
    FITS_OOSTan_FF5 = nan(size(RET,2),1);
    FITS_OOSTan_FF6 = nan(size(RET,2),1);

    %% Estimate

    oossuffix='';
    for n=1:size(RET,1)

        if oosflag
            for t=61:size(RET,2) % **** this changes across OOS specs
                if ~LOC(n,t) && Rlogic
                    continue
                end
                stind = 1;% **** this changes across OOS specs
                retnt = RET(n,stind:t-1)';
                ff1   = FF1(stind:t-1,:);
                ff3   = FF3(stind:t-1,:);
                ff4   = FF4(stind:t-1,:);
                ff5   = FF5(stind:t-1,:);
                ff6   = FF6(stind:t-1,:);
                ff1_mean = mean(ff1);%ff1_mean=mean(FF1);
                ff3_mean = mean(ff3);%ff3_mean=mean(FF3);
                ff4_mean = mean(ff4);%ff4_mean=mean(FF4);
                ff5_mean = mean(ff5);%ff5_mean=mean(FF5);
                ff6_mean = mean(ff6);%ff6_mean=mean(FF6);
                oossuffix = '_rec_60_60'; % what do you call this specification (that changed all the **** lines)
                if sum(~isnan(sum([retnt ff1 ff3 ff4 ff5 ff6],2)))<60% **** this changes across OOS specs
                    continue
                else
                    beta1 = regress(retnt,ff1);
                    beta3 = regress(retnt,ff3);
                    beta4 = regress(retnt,ff4);
                    beta5 = regress(retnt,ff5);
                    beta6 = regress(retnt,ff6);
                    FITS_FF1(n,t) = FF1(t,:)*beta1;
                    FITS_FF3(n,t) = FF3(t,:)*beta3;
                    FITS_FF4(n,t) = FF4(t,:)*beta4;
                    FITS_FF5(n,t) = FF5(t,:)*beta5;
                    FITS_FF6(n,t) = FF6(t,:)*beta6;
                    FITS_cond_FF1(n,t) = ff1_mean*beta1;
                    FITS_cond_FF3(n,t) = ff3_mean*beta3;
                    FITS_cond_FF4(n,t) = ff4_mean*beta4;
                    FITS_cond_FF5(n,t) = ff5_mean*beta5;
                    FITS_cond_FF6(n,t) = ff6_mean*beta6;

                    FITS_OOSTan_FF1(t)  = tanptfnext(ff6(:,1:1),FF6(t,1:1));
                    FITS_OOSTan_FF2(t)  = tanptfnext(ff6(:,1:2),FF6(t,1:2));
                    FITS_OOSTan_FF3(t)  = tanptfnext(ff6(:,1:3),FF6(t,1:3));
                    FITS_OOSTan_FF4(t)  = tanptfnext(ff6(:,1:4),FF6(t,1:4));
                    FITS_OOSTan_FF5(t)  = tanptfnext(ff6(:,1:5),FF6(t,1:5));
                    FITS_OOSTan_FF6(t)  = tanptfnext(ff6(:,1:6),FF6(t,1:6));
                end
            end
        else

            % Full-sample beta
            if sum(LOC(n,:),2)<12 && Rlogic
                continue
            end
            beta1 = regress(RET(n,:)',FF1);
            beta3 = regress(RET(n,:)',FF3);
            beta4 = regress(RET(n,:)',FF4);
            beta5 = regress(RET(n,:)',FF5);
            beta6 = regress(RET(n,:)',FF6);
            FITS_FF1(n,:) = FF1*beta1;
            FITS_FF3(n,:) = FF3*beta3;
            FITS_FF4(n,:) = FF4*beta4;
            FITS_FF5(n,:) = FF5*beta5;
            FITS_FF6(n,:) = FF6*beta6;
            FITS_cond_FF1(n,:) = FF1_mean*beta1;
            FITS_cond_FF3(n,:) = FF3_mean*beta3;
            FITS_cond_FF4(n,:) = FF4_mean*beta4;
            FITS_cond_FF5(n,:) = FF5_mean*beta5;
            FITS_cond_FF6(n,:) = FF6_mean*beta6;

        end

        if mod(n,100)==0;disp(['done with ' num2str(n)]);end
    end

    if Rlogic
    tmp                  = ~isnan(FITS_FF1+FITS_FF3+FITS_FF4+FITS_FF5+FITS_FF6);
    FITS_FFLOC           = false(size(LOC));
    FITS_FFLOC(:,loc2)   = tmp;
    else
    FITS_FFLOC           = ~isnan(FITS_FF1+FITS_FF3+FITS_FF4+FITS_FF5+FITS_FF6);
    end

    FITS_Factors = FF6;



    %% Save results

    if oosflag
        oosstring = 'OOS';
    else
        oosstring = '';
    end

    if ~Rlogic
    save(['../Data/Results_ObsFactReg' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    else
        save(['../Data/Results_ObsFactReg' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    end




%% SY
    case 'SY'
    %% Load data
    load(['../Data/' dataname]);

    if annualize
        ffsuffix = '_ANNRET';
    else
        ffsuffix = '';
    end
    sydata = load('../Data/M4');



    %% Start estimation
    disp(['Observable_Factor_Regressions starting at ' datestr(clock) ': data=' dataname]);



    %% FF models
    [date,loc1,loc2] = intersect(sydata.dates,date);
    SY = [sydata.MKTRF(loc1) sydata.SMB(loc1) sydata.MGMT(loc1) sydata.PERF(loc1)];

    switch QXorR
        case 'R'
            RET     = xret(:,loc2);
            Rlogic  = true;
        case 'Q'
            RET     = Q(:,loc2);
            Rlogic  = false;
        case 'X'
            RET     = X(:,loc2);
            Rlogic  = false;
    end

    % preallocate
    FITS_SY = nan(size(RET));

    FITS_cond_SY = nan(size(RET));


    SY_mean = repmat(mean(SY),size(RET,2),1);

    FITS_OOSTan_SY = nan(size(RET,2),1);


    %% Estimate

    oossuffix='';
    for n=1:size(RET,1)

        if oosflag
            for t=61:size(RET,2) % **** this changes across OOS specs
                if ~LOC(n,t) && Rlogic
                    continue
                end
                stind = 1;% **** this changes across OOS specs
                retnt = RET(n,stind:t-1)';
                sy   = SY(stind:t-1,:);
                sy_mean = mean(sy);%ff1_mean=mean(FF1);
                oossuffix = '_rec_60_60'; % what do you call this specification (that changed all the **** lines)
                if sum(~isnan(sum([retnt sy],2)))<60% **** this changes across OOS specs
                    continue
                else
                    beta1 = regress(retnt,sy);
                    FITS_SY(n,t) = SY(t,:)*beta1;
                    FITS_cond_SY(n,t) = sy_mean*beta1;

                    FITS_OOSTan_SY(t)  = tanptfnext(sy(:,1:1),SY(t,1:1));
                end
            end
        else

            % Full-sample beta
            if sum(LOC(n,:),2)<12 && Rlogic
                continue
            end
            beta1 = regress(RET(n,:)',SY);
            FITS_SY(n,:) = SY*beta1;
            FITS_cond_SY(n,:) = SY_mean*beta1;

        end

        if mod(n,100)==0;disp(['done with ' num2str(n)]);end
    end

    if Rlogic
    tmp                  = ~isnan(FITS_SY);
    FITS_SYLOC           = false(size(LOC));
    FITS_SYLOC(:,loc2)   = tmp;
    else
    FITS_SYLOC           = ~isnan(FITS_SY);
    end

    FITS_Factors = SY;



    %% Save results

    if oosflag
        oosstring = 'OOS';
    else
        oosstring = '';
    end

    if ~Rlogic
    save(['../Data/Results_ObsFactReg_SY_' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    else
        save(['../Data/Results_ObsFactReg_SY_' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    end
    
%% HXZ
    case 'HXZ'
    %% Load data
    load(['../Data/' dataname]);

    if annualize
        ffsuffix = '_ANNRET';
    else
        ffsuffix = '';
    end
    hxzdata = load('../Data/HXZ_q-Factors_monthly.mat');



    %% Start estimation
    disp(['Observable_Factor_Regressions starting at ' datestr(clock) ': data=' dataname]);



    %% FF models
    [date,loc1,loc2] = intersect(hxzdata.yrmo,date);
    HXZ = [hxzdata.Mkt_RF(loc1) hxzdata.ME(loc1) hxzdata.IA(loc1) hxzdata.ROE(loc1)];

    switch QXorR
        case 'R'
            RET     = xret(:,loc2);
            Rlogic  = true;
        case 'Q'
            RET     = Q(:,loc2);
            Rlogic  = false;
        case 'X'
            RET     = X(:,loc2);
            Rlogic  = false;
    end

    % preallocate
    FITS_HXZ = nan(size(RET));

    FITS_cond_HXZ = nan(size(RET));


    HXZ_mean = repmat(mean(HXZ),size(RET,2),1);

    FITS_OOSTan_HXZ = nan(size(RET,2),1);


    %% Estimate

    oossuffix='';
    for n=1:size(RET,1)

        if oosflag
            for t=61:size(RET,2) % **** this changes across OOS specs
                if ~LOC(n,t) && Rlogic
                    continue
                end
                stind = 1;% **** this changes across OOS specs
                retnt = RET(n,stind:t-1)';
                hxz   = HXZ(stind:t-1,:);
                hxz_mean = mean(hxz);%ff1_mean=mean(FF1);
                oossuffix = '_rec_60_60'; % what do you call this specification (that changed all the **** lines)
                if sum(~isnan(sum([retnt hxz],2)))<60% **** this changes across OOS specs
                    continue
                else
                    beta1 = regress(retnt,hxz);
                    FITS_HXZ(n,t) = HXZ(t,:)*beta1;
                    FITS_cond_HXZ(n,t) = hxz_mean*beta1;

                    FITS_OOSTan_HXZ(t)  = tanptfnext(hxz(:,1:1),HXZ(t,1:1));
                end
            end
        else

            % Full-sample beta
            if sum(LOC(n,:),2)<12 && Rlogic
                continue
            end
            beta1 = regress(RET(n,:)',HXZ);
            FITS_HXZ(n,:) = HXZ*beta1;
            FITS_cond_HXZ(n,:) = HXZ_mean*beta1;

        end

        if mod(n,100)==0;disp(['done with ' num2str(n)]);end
    end

    if Rlogic
    tmp                  = ~isnan(FITS_HXZ);
    FITS_HXZLOC          = false(size(LOC));
    FITS_HXZLOC(:,loc2)  = tmp;
    tmp                  = nan(size(LOC));
    tmp(:,loc2)          = FITS_HXZ;
    FITS_HXZ             = tmp;
    tmp                  = nan(size(LOC));
    tmp(:,loc2)          = FITS_cond_HXZ;
    FITS_cond_HXZ        = tmp;             
    else
    FITS_HXZLOC           = ~isnan(FITS_HXZ);
    end

    FITS_Factors = HXZ;



    %% Save results

    if oosflag
        oosstring = 'OOS';
    else
        oosstring = '';
    end

    if ~Rlogic
    save(['../Data/Results_ObsFactReg_HXZ_' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    else
        save(['../Data/Results_ObsFactReg_HXZ_' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    end
    
    
%% BS
    case 'BS'
    %% Load data
    load(['../Data/' dataname]);

    if annualize
        ffsuffix = '_ANNRET';
    else
        ffsuffix = '';
    end
    bsdata = load('../Data/BarillasShanken.mat');



    %% Start estimation
    disp(['Observable_Factor_Regressions starting at ' datestr(clock) ': data=' dataname]);



    %% FF models
    [date,loc1,loc2] = intersect(bsdata.dates,date);
    BS = [bsdata.Mkt_RF(loc1) bsdata.SMB(loc1) bsdata.UMD(loc1) bsdata.HMLm(loc1) bsdata.IA(loc1) bsdata.ROE(loc1)];

    switch QXorR
        case 'R'
            RET     = xret(:,loc2);
            Rlogic  = true;
        case 'Q'
            RET     = Q(:,loc2);
            Rlogic  = false;
        case 'X'
            RET     = X(:,loc2);
            Rlogic  = false;
    end

    % preallocate
    FITS_BS = nan(size(RET));

    FITS_cond_BS = nan(size(RET));


    BS_mean = repmat(mean(BS),size(RET,2),1);

    FITS_OOSTan_BS = nan(size(RET,2),1);


    %% Estimate

    oossuffix='';
    for n=1:size(RET,1)

        if oosflag
            for t=61:size(RET,2) % **** this changes across OOS specs
                if ~LOC(n,t) && Rlogic
                    continue
                end
                stind = 1;% **** this changes across OOS specs
                retnt = RET(n,stind:t-1)';
                bs   = BS(stind:t-1,:);
                bs_mean = mean(bs);%ff1_mean=mean(FF1);
                oossuffix = '_rec_60_60'; % what do you call this specification (that changed all the **** lines)
                if sum(~isnan(sum([retnt bs],2)))<60% **** this changes across OOS specs
                    continue
                else
                    beta1 = regress(retnt,bs);
                    FITS_BS(n,t) = BS(t,:)*beta1;
                    FITS_cond_BS(n,t) = bs_mean*beta1;

                    FITS_OOSTan_BS(t)  = tanptfnext(bs(:,1:1),BS(t,1:1));
                end
            end
        else

            % Full-sample beta
            if sum(LOC(n,:),2)<12 && Rlogic
                continue
            end
            beta1 = regress(RET(n,:)',BS);
            FITS_BS(n,:) = BS*beta1;
            FITS_cond_BS(n,:) = BS_mean*beta1;

        end

        if mod(n,100)==0;disp(['done with ' num2str(n)]);end
    end

    if Rlogic
    tmp                  = ~isnan(FITS_BS);
    FITS_BSLOC          = false(size(LOC));
    FITS_BSLOC(:,loc2)  = tmp;
    tmp                  = nan(size(LOC));
    tmp(:,loc2)          = FITS_BS;
    FITS_BS             = tmp;
    tmp                  = nan(size(LOC));
    tmp(:,loc2)          = FITS_cond_BS;
    FITS_cond_BS        = tmp;             
    else
    FITS_BSLOC           = ~isnan(FITS_BS);
    end

    FITS_Factors = BS;



    %% Save results

    if oosflag
        oosstring = 'OOS';
    else
        oosstring = '';
    end

    if ~Rlogic
    save(['../Data/Results_ObsFactReg_BS_' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    else
        save(['../Data/Results_ObsFactReg_BS_' QXorR oosstring '_' dataname oossuffix] ...
        , 'FITS*' , 'date' , 'RET')
    end
    

end