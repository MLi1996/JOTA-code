rng(2)

l = 100;
m = 100;
n = 200;


A = randn(m,n);

Phi = randn(l,n);
y = randn(l,1);
b = randn(m,1);



w = 0.1;

for i = 1 : m
    A(i,:) = A(i,:) / norm(A(i,:));
end

cvx_begin
    variable x(n,1);
    dual variable lambda;
    minimize (1/2/l * sum_square(Phi*x-y) + w / 2 * x'* x);
    subject to
        lambda: A * x <= b;
cvx_end

x_cvx = x;
lambda_cvx = lambda;

mul_SASC = 1;
eta_SASC = 2;

mul_SGD = 1.5;
eta_SGD = 2;

mul_SGDM = 1.5;
eta_SGDM = 2;
alpha_SGDM = 0.9;

mul_SVRG = 1;
eta_SVRG = 2;

mul_SVRG_screening = 1;
eta_SVRG_screening = 2;

Itrmax = 1e7 + 1;

sample = unidrnd(m,Itrmax,1); 

L = w + norm(Phi,2) ^ 2 / l;

mu = w;

Itrmax1 = Itrmax;
Itrmax2 = Itrmax;

sample2 = unidrnd(m,Itrmax2,1);

x0 = zeros(n,1);
xi = 5;
delta0 = 0.05;

SVRG_ratio = 3;
Catalyst_Budget = 1;


% Nested-Catalyst_SVRG

err_Catalyst_hist = zeros(Itrmax2,1);

ritr = 0;

x = x0;
mul = 2;
eta = 2;
N = SVRG_ratio;

T_B = Catalyst_Budget;
 
mtemp = m;

Itemp = 1:mtemp;

itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    xy = x;

    xold = x;

    inner_num = inner_num + 1;

    Ltemp = L + m * xi / 4 / delta;

    kappa = ( L - mu ) / m ;

    % kappa = 0;

    q = mu / ( mu + kappa );

    alpha = sqrt(q);

    beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );

    s = 3 / ( Ltemp + mu + kappa );

    M = ceil( mul * ( 16 * log(2*eta-1) + log( L / mu + mtemp * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + mtemp^2 * xi / 4 / mu / delta  ) );

    M = ceil( M / ( N + 1 ) / mtemp / T_B) * ( N + 1 ) * mtemp * T_B;

    sample_temp = unidrnd(mtemp,M,1);

    i = 0;
    t_b = 0;

    while i < M

        if mod(i,(N+1)*mtemp) == 0
            z = x;
            v1 = zeros(n,1);
            for j = 1 : l

                v1 = v1 + ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) / l;

            end
            
            ritr = ritr + l;

            for jj = 1 : mtemp
                
                j = Itemp(jj);
                v1 = v1 + xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

            end
                
            ritr = ritr + mtemp;

            i = i + mtemp;

        end

        ritr = ritr + 2;

        itr = itr + 1;
        i = i + 1;

        
        j = sample2(itr);
        j_cons = Itemp(sample_temp(i));
        

        
        gx1 = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) ;
        gz1 = ( Phi(j,:) )' * ( Phi(j,:) * z - y(j) ) ;
        
        gx2 = mtemp * xi * ( A(j_cons,:) )' * expd((A(j_cons,:)*x-b(j_cons))/delta);
        gz2 = mtemp * xi * ( A(j_cons,:) )' * expd((A(j_cons,:)*z-b(j_cons))/delta);

        gdiff = gx1 - gz1 + gx2 - gz2;

        x = ( x - s * gdiff - s * v1) / ( 1 + s * mu + s * kappa ) + xy * s * kappa / ( 1 + s * mu + s * kappa );

        err_Catalyst_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);

        if itr >= Itrmax2
            break
        end

        if mod(i,(N+1)*mtemp) == 0
            t_b = t_b + 1;
        end
        
        if t_b == T_B
            xy = x + beta * (x - xold);
            xold = x;
            t_b = 0;
        end
    end


       
    if itr >= Itrmax2
        break;
    end

%     Itemp_new = [];
% 
%     for i = 1 : mtemp
% 
%         j = Itemp(i);
% 
%         val = A(j,:)*x - b(j);
% 
%         if val >= -2 * sqrt(mtemp) * delta
%             Itemp_new = [Itemp_new,j];
%         end
%     end
% 
%     Itemp = Itemp_new;
%     mtemp = length(Itemp);

    delta = delta / eta;

end

inner_num_catalyst = inner_num;
ritr_Catalyst = ritr / 2;




% Nested-Catalyst_SVRG_Screening

err_Catalyst_Screening_hist = zeros(Itrmax2,1);

ritr = 0;

x = x0;
mul = 2;
eta = 2;
N = SVRG_ratio;

T_B = Catalyst_Budget;
 
mtemp = m;

Itemp = 1:mtemp;

itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    xy = x;

    xold = x;

    inner_num = inner_num + 1;

    Ltemp = L + m * xi / 4 / delta;

    kappa = ( L - mu ) / m ;

    % kappa = 0;

    q = mu / ( mu + kappa );

    alpha = sqrt(q);

    beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );

    s = 3 / ( Ltemp + mu + kappa );

    M = ceil( mul * ( 16 * log(2*eta-1) + log( L / mu + mtemp * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + mtemp^2 * xi / 4 / mu / delta  ) );

    M = ceil( M / ( N + 1 ) / mtemp / T_B) * ( N + 1 ) * mtemp * T_B;

    sample_temp = unidrnd(mtemp,M,1);

    i = 0;
    t_b = 0;

    while i < M

        if mod(i,(N+1)*mtemp) == 0
            z = x;
            v1 = zeros(n,1);
            for j = 1 : l

                v1 = v1 + ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) / l;

            end
            
            ritr = ritr + l;

            for jj = 1 : mtemp
                
                j = Itemp(jj);
                v1 = v1 + xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

            end
                
            ritr = ritr + mtemp;

            i = i + mtemp;

        end

        ritr = ritr + 2;

        itr = itr + 1;
        i = i + 1;

        
        j = sample2(itr);
        j_cons = Itemp(sample_temp(i));
        

        
        gx1 = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) ;
        gz1 = ( Phi(j,:) )' * ( Phi(j,:) * z - y(j) ) ;
        
        gx2 = mtemp * xi * ( A(j_cons,:) )' * expd((A(j_cons,:)*x-b(j_cons))/delta);
        gz2 = mtemp * xi * ( A(j_cons,:) )' * expd((A(j_cons,:)*z-b(j_cons))/delta);

        gdiff = gx1 - gz1 + gx2 - gz2;

        x = ( x - s * gdiff - s * v1) / ( 1 + s * mu + s * kappa ) + xy * s * kappa / ( 1 + s * mu + s * kappa );

        err_Catalyst_Screening_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);

        if itr >= Itrmax2
            break
        end

        if mod(i,(N+1)*mtemp) == 0
            t_b = t_b + 1;
        end
        
        if t_b == T_B
            xy = x + beta * (x - xold);
            xold = x;
            t_b = 0;
        end
    end


       
    if itr >= Itrmax2
        break;
    end

    Itemp_new = [];

    for i = 1 : mtemp

        j = Itemp(i);

        val = A(j,:)*x - b(j);

        if val >= -2 * sqrt(mtemp) * delta
            Itemp_new = [Itemp_new,j];
        end
    end

    Itemp = Itemp_new;
    mtemp = length(Itemp);

    delta = delta / eta;

end

inner_num_catalyst_screening = inner_num;
ritr_Catalyst_Screening = ritr / 2;

load('Exp1_SVRG_Screening.mat');
load('Exp1_SPDHG.mat');


% semilogy(    ritr_Catalyst_Screening / Itrmax2 * (1:100:Itrmax2), err_Catalyst_Screening_hist(1:100:Itrmax2), ...
%     1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax),...
%     (1+1/3)*(1:100:Itrmax),err_nested_SVRG_screening_hist(1:100:Itrmax));
% legend('Catalyst-Screening','SPDHG','Nested-SVRG-Screening');



semilogy(ritr_Catalyst / Itrmax2 * (1:100:Itrmax2), err_Catalyst_hist(1:100:Itrmax2),...
    ritr_Catalyst_Screening / Itrmax2 * (1:100:Itrmax2), err_Catalyst_Screening_hist(1:100:Itrmax2), ...
    1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax),...
    (1+1/3)*(1:100:Itrmax),err_nested_SVRG_screening_hist(1:100:Itrmax));
legend('Catalyst','Catalyst-Screening','SPDHG','Nested-SVRG-Screening');


% log_Err = zeros(m,1);
% 
% log_Err_Itemp = zeros(mtemp,1);
% 
% for k = 1 : m
% 
%     Err_temp = -(A(k,:)*x_cvx-b(k));
%     if Err_temp <= 1e-10
%         log_Err(k) = -10;
%     else
%         log_Err(k) = log10(Err_temp);
%     end
% 
% end
% 
% for i = 1 : mtemp
%     k = Itemp(i);
% 
%     Err_temp = -(A(k,:)*x_cvx-b(k));
%     if Err_temp <= 1e-10
%         log_Err_Itemp(i) = -10;
%     else
%         log_Err_Itemp(i) = log10(Err_temp);
%     end
% end
% 
% 
% edges = -10: 1: 2;
% h1 = histogram(log_Err,edges);
% hold on
% h2 = histogram(log_Err_Itemp,edges);
% legend('Original','Remaining')
% xlabel('lg(Slackness)');
% ylabel('Number of Constraints')