load("mushrooms_matrix.mat");
load("mushrooms_vector.mat");

m = 8124;
n = 111;

A = full(instance_matrix);

A = normalize(A);
A = [A(:,1:77),A(:,79:112)];

b = -ones(m,1);

for i = 1 : m
    nrm = norm(A(i,:),2);
    A(i,:) = A(i,:) / nrm;
    A(i,:) = - A(i,:) * label_vector(i);
    b(i) = b(i) / nrm;
end

rng(0)


xi = 1;
x0 = zeros(n,1);
delta0 = 0.005;

cvx_begin
    variable x(n,1);
    dual variable lambda;
    minimize 1 / 2 * (x'* x);
    subject to
        lambda: A * x <= b;
cvx_end

x_cvx = x;
lambda_cvx = lambda;

norm(lambda_cvx,'inf')


Itrmax = 1e7 + 1;

sample = unidrnd(m,Itrmax,1); 

L = 1;

mu = 1;

err_SASC_hist = zeros(Itrmax,1);
err_nested_SGD_hist = zeros(Itrmax,1);
err_nested_SGDM_hist = zeros(Itrmax,1);
err_static_SGDM_hist = zeros(Itrmax,1);


% SASC-SGD

x = x0;
mul = 1;
eta = 2;
itr = 0;

inner_num = 0;

s = 1 / 2; 

while inner_num >= 0

    inner_num = inner_num + 1;

    beta = 4 * s;

    M = ceil(mul * eta / mu / s );

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = (A(j,:))'*max(A(j,:)*x-b(j),0) / beta;

        x = ( x - s * g ) / ( 1 + s );

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
mul = 0.3;
eta = 2;
itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    inner_num = inner_num + 1;
    
    s = 4 / ( L + mu + m * xi / 4 / delta );

    M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

        x = ( x - s * g ) / ( 1 + s );

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

mul = 0.5;
eta = 4;
alpha = 0.5;
itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    inner_num = inner_num + 1;
    
    s = 2 / ( L + mu + m * xi / 4 / delta );

    M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + m * xi / 4 / mu / delta ) ) * sqrt( m * ( L / mu + m * xi / 4 / mu / delta ) ) );
    
    v = zeros(n,1);

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
        
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

s = 2 / ( L + mu + m * xi / 4 / delta );

for itr = 1 : Itrmax
    j = sample(itr);

    g = x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

    v = alpha * v - s * g;

    x = x + v ;

    err_static_SGDM_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
end
       

semilogy(1:100:Itrmax,err_SASC_hist(1:100:Itrmax),1:100:Itrmax,err_nested_SGD_hist(1:100:Itrmax),...
    1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax));
legend('SASC-SGD','Nested-SGD','Nested-SGDM','Static-SGDM');
xlabel('Iterations');
ylabel('Relative Error');