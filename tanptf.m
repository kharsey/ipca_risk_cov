function tp = tanptf(X)

% X is a matrix of excess returns, so fix rf=0

rf  = 0;
iota= ones(size(X,2),1);
S   = cov(X);
mu  = mean(X)';
tw  = S\(mu-rf*iota)/((iota'/S)*(mu-rf*iota));

tp  = X*tw;
