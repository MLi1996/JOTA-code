rng(3)

n = 500;
l = 100;
m = 100;

w_int = randn(n,1);

w_int = w_int / norm(w_int);

X10 = randn(l,n);
Y1 = zeros(l,1);

for i = 1:l
    X10(i,:) = X10(i,:) / norm(X10(i,:));
    Y1(i) = ( X10(i,:) * w_int > 0 );
end

X2 = randn(m,n);
Y2 = zeros(m,1);

for i = 1:m
    X2(i,:) = X2(i,:) / norm(X2(i,:));
    Y2(i) = ( X2(i,:) * w_int > 0 );
end

X11 = randn(l,n);


for i = 1:l
    X11(i,:) = X11(i,:) / norm(X11(i,:));
end

X1 = X10 + 0.1 * X11;


mu = 0.05;


A2 = zeros( m , n );

for i = 1 : m
    A2(i,:) = ( 2 * Y2(i) - 1 ) * X2(i,:);
end 


% Y1 = zeros(l,1);

cvx_begin
    variable w(n);
    dual variable lambda;
    minimize (- Y1' * X1 * w +sum(log_sum_exp([zeros(1,l); w'*X1'])) + mu / 2 * sum_square(w));
    subject to
        lambda: A2 * w >= 1;
cvx_end

w_cvx = w;
lambda_cvx = lambda;


w0 = zeros(n,1);
Itrmax = 1e5 + 1;
sample = unidrnd(l,Itrmax,1);

% SPDHG

err_SPDHG_hist = zeros(Itrmax,1);

% tau = 1 /  LA;
% 
% sigma = 1 / LA ;


Sigma = ones(m,2);

for i = 1 : m
    Sigma(i,2) = norm(X1(i,:),2);
end

gamma = 0.99;

tau = gamma / m / max(max(Sigma)) ;

sigma = gamma./ Sigma / 2 ;


w = w0;

psi = zeros(n,1);

for i = 1 : l
    psi = psi + Y1(i) * (X1(i,:))'; 
end

u = zeros(m,1);

v = zeros(l,1);

f = A2' * u;

f1 = f;

h = X1' * v;

h1 = h;

for i = 1 : Itrmax
    w = ( w - tau * f1 - tau * h1) / ( mu * tau + 1 ) +  psi * tau / ( 1 + tau * mu ) ;


    j = sample(i);
    utemp = min(0, u(j) + sigma(j,1) * (A2(j,:) * w - 1) );

    
    vtemp0 = v(j) + sigma(j,2) * X1(j,:) * w;

    func = @(val) sigma(j,2) * ( log(val) - log(1-val) ) + ( val - vtemp0 );

    options.Display = 'off';
    options.MaxIterations = 5;
    
    init = v(j);

    if init < 1e-10 || init > 1-1e-10
        init = 0.5;
    end

    vtemp = real(fsolve(func,init,options));
    
 
    theta = 1/sqrt(1+2*mu*tau);
    tau = theta * tau;
    sigma = sigma / theta;


    f = f + (A2(j,:))' * (utemp-u(j));

    f1 = f + (A2(j,:))' * (utemp-u(j)) * m * theta;
    
    h = h + (X1(j,:))' * (vtemp-v(j));

    h1 = h + (X1(j,:))' * (vtemp-v(j)) * l * theta;

    u(j) = utemp;
    
    v(j) = vtemp;

    err_SPDHG_hist(i) = norm(w-w_cvx,2)/norm(w_cvx,2);

end


semilogy(     1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax));
