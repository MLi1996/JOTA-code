% Multinomial Logistic Regression with Soft and Hard Classifications

rng(4)

num_C = 4;

n = 500;
l = 100;
m = 100;


W = randn(num_C,n);

for i = 1:num_C
    W(i,:) = W(i,:) / norm(W(i,:));
end

X10 = randn(n,l);
Z1 = zeros(1,l);
Y1 = zeros(num_C,l);

for i = 1:l
    X10(:,i) = X10(:,i) / norm(X10(:,i));
    [~, idxtemp] = max(W*X10(:,i));
    Z1(i) = idxtemp;
    Y1(idxtemp,i) = 1;
end

X2 = randn(n,m);
Z2 = zeros(1,m);
Y2 = zeros(num_C,m);

mytest = zeros(1,m);

for i = 1:m
    X2(:,i) = X2(:,i) / norm(X2(:,i));
    temp = W*X2(:,i);
    temp = sort(temp);
    mytest(i) = temp(4)-temp(3);
    [~, idxtemp] = max(W*X2(:,i));
    Z2(i) = idxtemp;
    Y2(idxtemp,i) = 1;
end

X11 = randn(n,l);


for i = 1:l
    X11(:,i) = X11(:,i) / norm(X11(:,i));
end

X1 = X10 + 0.1 * X11;


mu = 0.05;


A2 = zeros( (num_C-1) * m , num_C * n );

for i = 1 : m
    index = ( num_C - 1 ) * ( i - 1 );
    pos = Z2(i);
    pos_non = setdiff(1:num_C,pos); 
    for j = 1 : num_C - 1
        A2(index + j, ( (pos-1) * n + 1 ) : (pos * n) ) = (X2(:,i))';
        posn = pos_non(j);
        A2(index + j, ( (posn-1) * n + 1 ) : (posn * n) ) = -(X2(:,i))';
    end
end


cvx_begin
    variable w(n,num_C);
    dual variable lambda;
    minimize (-(vec(w*Y1))'*(vec(X1)) + sum(log_sum_exp((w')*X1)) + mu / 2 * sum_square(vec(w)));
    subject to
        lambda: A2 * vec(w) >= 1;
cvx_end

w_cvx = w;
lambda_cvx = lambda;


cvx_begin
    variable w(n,num_C);
    dual variable lambda;
    minimize (-(vec(w*Y1))'*(vec(X1)) + sum(log_sum_exp((w')*X1)) + mu / 2 * sum_square(vec(w)));
cvx_end

w_cvx2 = w;

norm(w_cvx-w_cvx2) / norm(w_cvx)


norm(lambda_cvx,'inf')



mul_SASC = 1;
eta_SASC = 2;

mul_SGD = 1;
eta_SGD = 2;

mul_SGDM = 1;
eta_SGDM = 2;
alpha_SGDM = 0.5;

%%%%%%


Itrmax = 1e6 + 1;

sample = unidrnd(m,Itrmax,1); 
% sample2 = unidrnd(m,Itrmax,1);

S = kron((eye(num_C)-ones(num_C,num_C)/num_C), X1*X1') / 2;

L = max(eig(S)) + mu;

mu = mu;

err_SASC_hist = zeros(Itrmax,1);
err_nested_SGD_hist = zeros(Itrmax,1);
err_nested_SGDM_hist = zeros(Itrmax,1);
err_static_SGDM_hist = zeros(Itrmax,1);

xi = 1;
delta0 = 1;


w0 = zeros(n,num_C);

% SASC-SGD

w = w0;
mul = mul_SASC;
eta = eta_SASC;
itr = 0;

inner_num = 0;

s = 1 / 2 / L; 

while inner_num >= 0

    inner_num = inner_num + 1;

    beta = 4 * s;

    M = ceil(mul * eta / mu / s );

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = mu * w;

        g(:,Z1(j)) = g(:,Z1(j)) - l * X1(:,j);

        myexp1 = exp(w'*X1(:,j));

        for k = 1 : num_C
            g(:,k) = g(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j); 
        end
        
        pos = Z2(j);

        pos_non = setdiff(1:num_C,pos);


        for k = 1 : num_C - 1
            val = max((X2(:,j))'*(w(:,pos_non(k))-w(:,pos))+1,0) / beta;
            g(:,pos_non(k)) = g(:,pos_non(k)) + val * X2(:,j);
            g(:,pos) = g(:,pos) - val * X2(:,j);
        end

        % g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + (A(j,:))'*max(A(j,:)*x-b(j),0) / beta;

        w = w - s * g;

        err_SASC_hist(itr) = norm(w-w_cvx,'fro')/norm(w_cvx,'fro');

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

w = w0;
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
        
        g = mu * w;

        g(:,Z1(j)) = g(:,Z1(j)) - l * X1(:,j);

        myexp1 = exp(w'*X1(:,j));

        for k = 1 : num_C
            g(:,k) = g(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j); 
        end
        
        pos = Z2(j);

        pos_non = setdiff(1:num_C,pos);

        for k = 1 : num_C - 1
            val = expd( ((X2(:,j))' * ( w(:,pos_non(k)) - w(:,pos) )+1) / delta ) * m * xi;
            g(:,pos_non(k)) = g(:,pos_non(k)) + val * X2(:,j);
            g(:,pos) = g(:,pos) - val * X2(:,j);
        end

        w = w - s * g;

        err_nested_SGD_hist(itr) = norm(w-w_cvx,'fro')/norm(w_cvx,'fro');

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

w = w0;

v = zeros(n,num_C);

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
    
    v = zeros(n,num_C);

    for i = 1 : M
        itr = itr + 1;
        j = sample(itr);
        
        g = mu * w;

        g(:,Z1(j)) = g(:,Z1(j)) - l * X1(:,j);

        myexp1 = exp(w'*X1(:,j));

        for k = 1 : num_C
            g(:,k) = g(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j);
        end

        pos = Z2(j);

        pos_non = setdiff(1:num_C,pos);

        for k = 1 : num_C - 1
            val = expd( ((X2(:,j))' * ( w(:,pos_non(k)) - w(:,pos) )+1) / delta ) * m * xi;
            g(:,pos_non(k)) = g(:,pos_non(k)) + val * X2(:,j);
            g(:,pos) = g(:,pos) - val * X2(:,j);
        end

        v = alpha * v - s * g;

        w = w + v;

        err_nested_SGDM_hist(itr) = norm(w-w_cvx,'fro')/norm(w_cvx,'fro');

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

w = w0;

v = zeros(n,num_C);

delta = sqrt(delta_nested_SGDM *delta0);

s = 1 / ( L + mu + m * xi / 4 / delta );

for itr = 1 : Itrmax
    j = sample(itr);

    g = mu * w;

    g(:,Z1(j)) = g(:,Z1(j)) - l * X1(:,j);

    myexp1 = exp(w'*X1(:,j));

    for k = 1 : num_C
        g(:,k) = g(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j);
    end

    pos = Z2(j);

    pos_non = setdiff(1:num_C,pos);

    for k = 1 : num_C - 1
        val = expd( ((X2(:,j))' * ( w(:,pos_non(k)) - w(:,pos) )+1) / delta ) * m * xi;
        g(:,pos_non(k)) = g(:,pos_non(k)) + val * X2(:,j);
        g(:,pos) = g(:,pos) - val * X2(:,j);
    end


    v = alpha * v - s * g;

    w = w + v;

    err_static_SGDM_hist(itr) = norm(w-w_cvx,'fro')/norm(w_cvx,'fro');
end




semilogy(1:100:Itrmax,err_SASC_hist(1:100:Itrmax),1:100:Itrmax,err_nested_SGD_hist(1:100:Itrmax),...
    1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax));
legend('SASC-SGD','Nested-SGD','Nested-SGDM','Static-SGDM');
ylabel('Relative Error');
xlabel('Iterations');