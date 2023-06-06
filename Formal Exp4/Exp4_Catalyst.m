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


cvx_begin
    variable w(n);
    dual variable lambda;
    minimize (- Y1' * X1 * w +sum(log_sum_exp([zeros(1,l); w'*X1'])) + mu / 2 * sum_square(w));
    subject to
        lambda: A2 * w >= 1;
cvx_end

w_cvx = w;
lambda_cvx = lambda;


norm(lambda_cvx,'inf')

%%%%%%


Itrmax = 1e6 + 1;

sample = unidrnd(m,Itrmax,1); 

num_C = 2;

S = kron((eye(num_C)-ones(num_C,num_C)/num_C), X1*X1') / 2;

L = max(eig(S)) + mu;

mu = mu;

xi = 1;
delta0 = 1;

SVRG_ratio = 3;
Catalyst_Budget = 1;

Itrmax1 = Itrmax;
Itrmax2 = Itrmax;

sample2 = unidrnd(m,Itrmax2,1);



err_Catalyst_hist = zeros(Itrmax2,1);

% Nested-Catalyst_SVRG

ritr = 0;

w0 = zeros(n,1);

w = w0;
mul = 1.5;
eta = 2;
N = SVRG_ratio;

T_B = Catalyst_Budget;
 
mtemp = m;

Itemp = 1:mtemp;

itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    wy = w;

    wold = w;

    inner_num = inner_num + 1;

    Ltemp = L + m * xi / 4 / delta;

    kappa = ( L - mu ) / m ;

    q = mu / ( mu + kappa );

    alpha = sqrt(q);

    beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );

    s = 3 / ( Ltemp + mu + kappa );

    M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + mtemp * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + mtemp^2 * xi / 4 / mu / delta  ) );

    M = ceil( M / ( N + 1 ) / mtemp / T_B) * ( N + 1 ) * mtemp * T_B;

    sample_temp = unidrnd(mtemp,M,1);

    i = 0;
    t_b = 0;

    while i < M

        if mod(i,(N+1)*mtemp) == 0
            y = w;
            v1 = zeros(n,1);
            for j = 1 : m

                g = - (Y1(j)*X1(j,:))';
                g = g + expd(X1(j,:)*w) * (X1(j,:))';
                v1 = v1 + g;

            end
            ritr = ritr + m;

            for jj = 1 : mtemp
                
                j = Itemp(jj);
                val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
                v1 = v1 - val * (A2(j,:))';

            end
                
            ritr = ritr + mtemp;
            i = i + mtemp;

        end

        ritr = ritr + 2;

        itr = itr + 1;
        i = i + 1;

        
        j = sample2(itr);
        j_cons = Itemp(sample_temp(i));

        gx = - l * (Y1(j)*X1(j,:))';

        gx = gx + l * expd(X1(j,:)*w) * (X1(j,:))';

        val = expd( (-A2(j_cons,:) * w + 1) / delta ) * mtemp * xi;

        gx = gx - val * (A2(j_cons,:))';

        gy = - l * (Y1(j)*X1(j,:))';

        gy = gy + l * expd(X1(j,:)*y) * (X1(j,:))';

        val = expd( (-A2(j_cons,:) * y + 1) / delta ) * mtemp * xi;

        gy = gy - val * (A2(j_cons,:))';

        gdiff = gx - gy;

        w = ( w - s * gdiff - s * v1) / ( 1 + s * mu + s * kappa ) + wy * s * kappa / ( 1 + s * mu + s * kappa );

        err_Catalyst_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);

        if itr >= Itrmax2
            break
        end

        if mod(i,(N+1)*mtemp) == 0
            t_b = t_b + 1;
        end
        
        if t_b == T_B
            wy = w + beta * (w - wold);
            wold = w;
            t_b = 0;
        end
    end


       
    if itr >= Itrmax2
        break;
    end
% 
%     Itemp_new = [];
% 
%     for i = 1 : mtemp
% 
%         j = Itemp(i);
% 
%         val = -A2(j,:) * w + 1;
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

ritr_Catalyst = ritr / 2;

inner_num_catalyst = inner_num;




err_Catalyst_Screening_hist = zeros(Itrmax2,1);

% Nested-Catalyst_SVRG_Screening

ritr = 0;

w0 = zeros(n,1);

w = w0;
mul = 1.5;
eta = 2;
N = SVRG_ratio;

T_B = Catalyst_Budget;
 
mtemp = m;

Itemp = 1:mtemp;

itr = 0;

inner_num = 0;

delta = delta0;

while inner_num >= 0

    wy = w;

    wold = w;

    inner_num = inner_num + 1;

    Ltemp = L + m * xi / 4 / delta;

    kappa = ( L - mu ) / m ;

    q = mu / ( mu + kappa );

    alpha = sqrt(q);

    beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );

    s = 3 / ( Ltemp + mu + kappa );

    M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + mtemp * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + mtemp^2 * xi / 4 / mu / delta  ) );

    M = ceil( M / ( N + 1 ) / mtemp / T_B) * ( N + 1 ) * mtemp * T_B;

    sample_temp = unidrnd(mtemp,M,1);

    i = 0;
    t_b = 0;

    while i < M

        if mod(i,(N+1)*mtemp) == 0
            y = w;
            v1 = zeros(n,1);
            for j = 1 : m

                g = - (Y1(j)*X1(j,:))';
                g = g + expd(X1(j,:)*w) * (X1(j,:))';
                v1 = v1 + g;

            end
            
            ritr = ritr + m;

            for jj = 1 : mtemp
                
                j = Itemp(jj);
                val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
                v1 = v1 - val * (A2(j,:))';

            end
                
            ritr = ritr + mtemp;

            i = i + mtemp;

        end

        ritr = ritr + 2;

        itr = itr + 1;
        i = i + 1;

        
        j = sample2(itr);
        j_cons = Itemp(sample_temp(i));

        gx = - l * (Y1(j)*X1(j,:))';

        gx = gx + l * expd(X1(j,:)*w) * (X1(j,:))';

        val = expd( (-A2(j_cons,:) * w + 1) / delta ) * mtemp * xi;

        gx = gx - val * (A2(j_cons,:))';

        gy = - l * (Y1(j)*X1(j,:))';

        gy = gy + l * expd(X1(j,:)*y) * (X1(j,:))';

        val = expd( (-A2(j_cons,:) * y + 1) / delta ) * mtemp * xi;

        gy = gy - val * (A2(j_cons,:))';

        gdiff = gx - gy;

        w = ( w - s * gdiff - s * v1) / ( 1 + s * mu + s * kappa ) + wy * s * kappa / ( 1 + s * mu + s * kappa );

        err_Catalyst_Screening_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);

        if itr >= Itrmax2
            break
        end

        if mod(i,(N+1)*mtemp) == 0
            t_b = t_b + 1;
        end
        
        if t_b == T_B
            wy = w + beta * (w - wold);
            wold = w;
            t_b = 0;
        end
    end


       
    if itr >= Itrmax2
        break;
    end

    Itemp_new = [];

    for i = 1 : mtemp

        j = Itemp(i);

        val = -A2(j,:) * w + 1;

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

load('Exp4_SPDHG.mat');
load('Exp4_SVRG_Screening.mat');

semilogy(ritr_Catalyst / Itrmax2 * (1:100:Itrmax2), err_Catalyst_hist(1:100:Itrmax2),...
    ritr_Catalyst_Screening / Itrmax2 * (1:100:Itrmax2), err_Catalyst_Screening_hist(1:100:Itrmax2),1.33333*(1:100:Itrmax2),err_nested_SVRG_Screening_hist(1:100:Itrmax2),...
     1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax));
legend('Nested-Catalyst','Nested-Catalyst-Screening','Nested-SVRG-Screening','SPDHG');
xlabel('Number of Stochastic Evaluations');
ylabel('Relative Error');



% log_Err = zeros(m,1);
% 
% log_Err_Itemp = zeros(mtemp,1);
% 
% for k = 1 : m
% 
%     Err_temp = (A2(k,:)*w_cvx-1);
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
%     Err_temp = (A2(k,:)*w_cvx-1);
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