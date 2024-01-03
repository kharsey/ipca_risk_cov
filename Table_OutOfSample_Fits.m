
clear


dirname = '../Data/';

name1 = 'Results_GB_outofsample_IPCADATA_FNW36_RNKDMN_CON_K'
nameFFR = 'Results_ObsFactRegROOS_Results_GB_IPCADATA_FNW36_RNKDMN_CON_K1_rec_60_60'
nameFFX = 'Results_ObsFactRegXOOS_Results_GB_IPCADATA_FNW36_RNKDMN_CON_K1_rec_60_60'

Krange = [1:6];

%% IPCA

si=120;
report=nan(30);% make it bigger than it needs to be -- only specific entries are displayed

for K=Krange
    try 
        load([dirname name1 num2str(K)]);
        ffr     = load([dirname nameFFR]);
        ffx     = load([dirname nameFFX]);
    catch
        continue
    end
    
    BIGLOC              = LOC & ffr.FITS_FFLOC;
    BIGLOC(:,1:si-1)    = false;
    BIGLOC = LOC;
    BIGLOC(:,1:si-1)    = false;
    
    try
    report(1,K) = ...
        1 - mySOS(xret(BIGLOC(:))-OOSRFITS_GB(BIGLOC(:)))/mySOS(xret(BIGLOC(:)));
    report(2,K) = ...
        1 - mySOS(xret(BIGLOC(:))-OOSRFITS_pred_GB(BIGLOC(:)))/mySOS(xret(BIGLOC(:)));
    report(3,K) = ...
        1 - mySOS(X(:,si:end)-OOSXFITS_GB(:,si:end))/mySOS(X(:,si:end));
    report(4,K) = ...
        1 - mySOS(X(:,si:end)-OOSXFITS_pred_GB(:,si:end))/mySOS(X(:,si:end));
    
    report(5,K) = ...
        1 - mySOS(xret(BIGLOC(:))-ffr.(['FITS_FF' num2str(K)])(BIGLOC(:)))/mySOS(xret(BIGLOC(:)));
    report(6,K) = ...
        1 - mySOS(xret(BIGLOC(:))-ffr.(['FITS_cond_FF' num2str(K)])(BIGLOC(:)))/mySOS(xret(BIGLOC(:)));
%     report(7,K) = ...
%         1 - mySOS(X(:,si:end)-ffx.(['FITS_FF' num2str(K)])(:,si:end))/mySOS(X(:,si:end));
%     report(8,K) = ...
%         1 - mySOS(X(:,si:end)-ffx.(['FITS_cond_FF' num2str(K)])(:,si:end))/mySOS(X(:,si:end));
    catch
    end
    

    tmp = OOSRealFact(:,si:end);
    report(9,K) = ...
        sqrt(12)*mean(tmp(end,:))/std(tmp(end,:));
    report(10,K) = ...
        sqrt(12)*mean(OOSRealTan(si:end))/std(OOSRealTan(si:end));
    report(11,K) = ...
        sqrt(12)*mean(ffr.FITS_Factors(si:end,K))/std(ffr.FITS_Factors(si:end,K));
    tmp = tanptf(ffr.FITS_Factors(:,1:K));
    report(12,K) = ...  
        sqrt(12)*mean(tmp(si:end))/std(tmp(si:end));
        
    
    
end

reportrows = {  'R_R2_total_GB (IPCA)        ' ...
                'R_R2_pred_GB (IPCA)         ' ...
                'X_R2_total_GB (IPCA)     ' ...
                'X_R2_pred_GB (IPCA)      ' ...
                'R_R2_total_GB (FF)         ' ...
                'R_R2_pred_GB (FF)          ' ...
                'X_R2_total_GB (FF)      ' ...
                'X_R2_pred_GB (FF)       ' ...
                'Univar Ptf (IPCA)           ' ...
                'Tangency Ptf (IPCA)         ' ...
                'Univar Ptf (FF)            ' ...
                'Tangency Ptf (FF)          '};

disp(' ')
disp([          'Ks :                   ' num2str(Krange,'%0.0f      ')])

for j=[1:length(reportrows)]
    rstring = [reportrows{j} ' '];
    for jj=Krange
        tmp = num2str(100*report(j,jj),3);
        rstring = [rstring ' & ' tmp];
    end
    rstring = [rstring ' \\'];
    disp(rstring)
end

disp(' ')
disp(char(reportrows'))
report(1:length(reportrows),Krange)




