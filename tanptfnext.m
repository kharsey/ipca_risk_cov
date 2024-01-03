function [tpnext,tw] = tanptfnext(X,Xnext,varargin)

% X     is a matrix of excess returns, so fix rf=0
% Xnext is a vector of future excess returns to apply weights to

rf      = 0;
iota    = ones(size(X,2),1);
S       = cov(X);
mu      = mean(X)';
tw      = S\(mu-rf*iota)/((iota'/S)*(mu-rf*iota));
if tw'*(mu-rf)<0, tw = -tw; end        

if nargin>=3
    targetvol   = varargin{1};
    scale = targetvol/sqrt(tw'*S*tw);
else
    scale = 1;
end

tw = tw*scale;

tpnext  = Xnext*tw;

end