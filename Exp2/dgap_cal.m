function [val] = dgap_cal(H,c,A,b,G,x,d,xi,delta)
% Calculation of the Duality Gap

val = 1 / 2 * x' * ( H * x )+ c' * x ; 

[m,~] = size(A);
lambda = zeros(m,1);

for i = 1 : m
    s = (A(i,:) * x - b(i))/delta;
    val = val + xi * delta * log( 1 + exp( s ) );
    if isnan(exp(s)) || exp(s) > 1e10
        lambda(i) = xi;
    else
        lambda(i) = xi * exp(s) / ( 1 + exp(s) );
    end
end

val = val + b' * lambda + 1 / 2 * c' * d + d' * A' * lambda + 1 / 2 * lambda' * (G *lambda);

end