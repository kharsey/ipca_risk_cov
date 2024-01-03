% Arb portfolios
clear
for j=1:6
    load(['../Data/Results_GBGA_outofsample_IPCADATA_FNW36_RNKDMN_CON_K' num2str(j) '.mat'])
    loc = floor(length(OOSARBPTF)/2):length(OOSARBPTF);
    disp([j nanmean(OOSARBPTF(loc))/nanstd(OOSARBPTF(loc))*sqrt(12)])
end