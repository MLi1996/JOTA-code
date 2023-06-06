function [U] = newtransform(x,K,k)
% Find an orthogonal transformation, s.t.
% UK = K, and Ux in span (K, k)

K2 = K'*K;

xres = x - K * ( K2 \ ( K' * x ) );

kres = k - K * ( K2 \ ( K' * k ) );

if norm(xres,2) < 1e-16
    U = eye(length(x));
    return
end

xres = xres / norm(xres,2);

kres = kres / norm(kres,2);

w = xres - kres;

w = w / norm(w,2);

U = eye(length(x)) - 2 * (w * w');

end

