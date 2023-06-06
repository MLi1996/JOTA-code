function [val] = expd(t)
%

s = exp(t);

if isnan(s) || s > 1e10
    val = 1;
else
    val = s / ( 1 + s );
end



end