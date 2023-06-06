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

H = eye(n);

G = A * (H \ A');

c = zeros(n,1);

d = c;


Itrmax = 1e7 + 1;

sample = unidrnd(m,Itrmax,1); 

L = 1;

mu = 1;

dgap_nested_SGD_hist = zeros(Itrmax,1);
dgap_nested_SGDM_hist = zeros(Itrmax,1);
dgap_static_SGDM_hist = zeros(Itrmax,1);

dgap_ns_min = 1e9;
dgap_nm_min = 1e9;
dgap_sm_min = 1e9;


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

        if mod(itr,1000) == 1
            dtemp = dgap_cal(H,c,A,b,G,x,d,xi,delta);
            if dgap_ns_min > dtemp
                dgap_ns_min = dtemp;
            end
            dgap_nested_SGD_hist(itr) = dgap_ns_min;
        end

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

        if mod(itr,1000) == 1
            dtemp = dgap_cal(H,c,A,b,G,x,d,xi,delta);
            if dgap_nm_min > dtemp
                dgap_nm_min = dtemp;
            end
            dgap_nested_SGDM_hist(itr) = dgap_nm_min;
        end  

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

    if mod(itr,1000) == 1
        dtemp = dgap_cal(H,c,A,b,G,x,d,xi,delta);
        if dgap_sm_min > dtemp
            dgap_sm_min = dtemp;
        end
        dgap_static_SGDM_hist(itr) = dgap_sm_min;
    end
end
       

semilogy(1:1000:Itrmax,dgap_nested_SGD_hist(1:1000:Itrmax),...
    1:1000:Itrmax,dgap_nested_SGDM_hist(1:1000:Itrmax),1:1000:Itrmax,dgap_static_SGDM_hist(1:1000:Itrmax));
legend('Nested-SGD','Nested-SGDM','Static-SGDM');
xlabel('Iterations');
ylabel('Minimal Duality Gap');