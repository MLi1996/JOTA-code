function [g] = newgrad(x, H, xi, delta, A, b)
% Calculate the gradient

g = H * x;

m = length(b);

for i = 1 : m
    g = g + xi * (exp((A(i,:)*x-b(i))/delta) - 1) / (exp((A(i,:)*x-b(i))/delta) + 1) * (A(i,:))';
end

end

