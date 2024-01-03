function ret =  mySOS(p)
    p  = p(:);
    po = ~isnan(p);
    ret = p(po(:))'*p(po(:));
%     ret = trace(p * p');
end