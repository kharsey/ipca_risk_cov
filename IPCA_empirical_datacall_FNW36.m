%% this data has characteristics lagged (ARE LAGGED)
load('../Data/characteristics_data_feb2017.mat');

chars = char;
date = yrmo;

% make lme (lagged me) into log me
loc = strmatch('lme',charnames,'exact');
chars(:,:,loc) = log(chars(:,:,loc));

% make at into log at
loc = strmatch('at',charnames,'exact');
chars(:,:,loc) = log(chars(:,:,loc));

%% make xret
ffd = load('../Data/FF3F_RF_192607_201606.mat');
[~,loc1,loc2] = intersect(date,ffd.yrmolistFF);
if length(find(loc1))~= length(date)
    error('The provided RF series does not completely overlap the stock/chars data sample')
end
xret = bsxfun( @minus , ret , ffd.RF(loc2)'/100 );
RF   = ffd.RF(loc2)'/100;

%% PUBLICATION ORDER (for FNW36) 
puborder = ...
{   'a2me'                      1988;
    'at'                        2015;
    'ato'                       2008;
    'beme'                      1985;
    'beta'                      1973; % ref 2014 paper, but I go back to Black Jensen Scholes
    'c'                         2012;
    'cto'                       1996;
    'd2a'                       2016;
    'dpi2a'                     2008;
    'e2p'                       1983;
    'fc2y'                      2016;
    'free_cf'                   2011;
    'idio_vol'                  2006;
    'investment'                2008;
    'lev'                       2015;
    'lme'                       1992;% this is mistaken, should be banz 1985
    'lturnover'                 1998;
    'noa'                       2004;
    'oa'                        1996;
    'ol'                        2011;
    'pcm'                       2016;
    'pm'                        2008;
    'prof'                      2015;
    'q'                         1985; % don't state a date, so I put it in the same year as beme
    'rel_to_high_price'         2004;
    'rna'                       2008;
    'roa'                       2010;
    'roe'                       1996;
    'cum_return_12_2'           1996;
    'cum_return_12_7'           2012;
    'cum_return_1_0'            1990;
    'cum_return_36_13'          1985;
    's2p'                       2015;
    'sga2m'                     2015; % is misnamed "sga2m" instead of "sga2s", and don't state a date, so I put with Lewellen's 2015
    'spread_mean'               2014;
    'suv'                       2009};
[~,pubordersort] = sort([puborder{:,2}]);
puborder = puborder(pubordersort,:);
[~,~,L] = size(chars);
tmp = nan(L,1);
for l=1:L
    tmp(l) = find( strcmp(puborder{l,1},charnames) );
end
chars = chars(:,:,tmp);
charnames = charnames(tmp);
clear tmp

clear char yrmo loc ffd ret
