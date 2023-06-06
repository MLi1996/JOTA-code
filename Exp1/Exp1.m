rng(2)

l = 100;
m = 100;
n = 100;


A = randn(m,n);

Phi = randn(l,n);
y = randn(l,1);
b = abs(randn(m,1));





w = 0.1;

xi = 1;

x0 = zeros(n,1);

for i = 1 : m
    A(i,:) = A(i,:) / norm(A(i,:));
end

delta0 = 0.05;

cvx_begin
    variable x(n,1);
    dual variable lambda;
    minimize (1/2/l * sum_square(Phi*x-y) + w / 2 * x'* x);
    subject to
        lambda: A * x <= b;
cvx_end

x_cvx = x;
lambda_cvx = lambda;

norm(lambda_cvx,'inf')

mul_SASC = 0.36;
eta_SASC = 2;

mul_SGD = 0.6;
eta_SGD = 4;

mul_SGDM = 1;
eta_SGDM = 2;
alpha_SGDM = 0.9;

Itrmax = 1e7 + 1;

sample = unidrnd(m,Itrmax,1); 

L = w + norm(Phi,2) ^ 2 / l;

mu = w;

err_SASC_hist = zeros(Itrmax,1);
err_nested_SGD_hist = zeros(Itrmax,1);
err_nested_SGDM_hist = zeros(Itrmax,1);
err_static_SGDM_hist = zeros(Itrmax,1);


% SASC-SGD

x = x0;
mul = mul_SASC;
eta = eta_SASC;
itr = 0;

inner_num = 0;

s = 1 / 6 / L; 

while inner_num >= 0

    inner_num = inner_num + 1;

    beta = 4 * s;

    M = ceil(mul * eta / mu / s );

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + (A(j,:))'*max(A(j,:)*x-b(j),0) / beta;

        x = x - s * g;

        err_SASC_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);

        if itr >= Itrmax
            break
        end

    end

    if itr >= Itrmax
        break
    end

    s = s / eta;

end

inner_num_SASC = inner_num;

% Nested-SGD

x = x0;
mul = mul_SGD;
eta = eta_SGD;
itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    inner_num = inner_num + 1;
    
    s = 1 / ( L + mu + m * xi / 4 / delta );

    M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

        x = x - s * g;

        err_nested_SGD_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);

        if itr >= Itrmax
            break
        end

    end

    if itr >= Itrmax
        break
    end

    delta = delta / eta;

end

inner_num_SGD = inner_num;

% Nested-SGDM

x = x0;

v = zeros(n,1);

mul = mul_SGDM;
eta = eta_SGDM;
alpha = alpha_SGDM;
itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    inner_num = inner_num + 1;
    
    s = 1 / ( L + mu + m * xi / 4 / delta );

    M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + m * xi / 4 / mu / delta ) ) * sqrt( m * ( L / mu + m * xi / 4 / mu / delta ) ) );
    
    v = zeros(n,1);

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
        
        v = alpha * v - s * g;

        x = x + v;

        err_nested_SGDM_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);

        if itr >= Itrmax
            break
        end

    end

    if itr >= Itrmax
        break
    end

    delta = delta / eta;

end

inner_num_SGDM = inner_num;

delta_nested_SGDM = delta;

% Static-SGDM

x = x0;

v = zeros(n,1);

delta = sqrt(delta_nested_SGDM *delta0);

s = 1 / ( L + mu + m * xi / 4 / delta );

for itr = 1 : Itrmax
    j = sample(itr);

    g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

    v = alpha * v - s * g;

    x = x + v;

    err_static_SGDM_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
end


% SPDHG

err_SPDHG_hist = zeros(Itrmax,1);

LA = sqrt(m +  l * (L - mu));

tau = 1 /  LA;

sigma = 1 / LA ;

x = x0;

% u = xi*rand(m,1);

u = 0.5*rand(m,1);
v = zeros(l,1);
f = A' * u;
h = Phi' * v;

for i = 1 : Itrmax
    x = ( x - tau * f - tau * h) / ( w * tau + 1 );


    j = sample(i);
    utemp = u(j) + sigma * A(j,:) * x;
    utemp2 = max(0,utemp-sigma*b(j));
    
    vtemp = (v(j)+sigma*(Phi(j,:)*x-y(j))) / (1+sigma*l);

    theta = 1/sqrt(1+2*mu*tau);
    tau = theta * tau;
    sigma = sigma / theta;


    f = f + (A(j,:))' * (utemp2-u(j));

    f = f + (A(j,:))' * (utemp2-u(j)) / m * theta;
    
    h = h + (Phi(j,:))' * (vtemp-v(j));

    h = h + (Phi(j,:))' * (vtemp-v(j)) / l * theta;

    u(j) = utemp2;
    
    v(j) = vtemp;

    err_SPDHG_hist(i) = norm(x-x_cvx,2)/norm(x_cvx,2); 
end



semilogy(1:100:Itrmax,err_SASC_hist(1:100:Itrmax),1:100:Itrmax,err_nested_SGD_hist(1:100:Itrmax),...
    1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax),...
    1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax));
legend('SASC-SGD','Nested-SGD','Nested-SGDM','Static-SGDM','SPDHG');
