% rng(2)
% 
% l = 100;
% m = 100;
% n = 200;
% 
% 
% A = randn(m,n);
% 
% Phi = randn(l,n);
% y = randn(l,1);
% b = randn(m,1);
% 
% 
% 
% w = 0.1;
% 
% xi = 5;
% 
% x0 = zeros(n,1);
% 
% for i = 1 : m
%     A(i,:) = A(i,:) / norm(A(i,:));
% end
% 
% delta0 = 0.05;
% 
% cvx_begin
%     variable x(n,1);
%     dual variable lambda;
%     minimize (1/2/l * sum_square(Phi*x-y) + w / 2 * x'* x);
%     subject to
%         lambda: A * x <= b;
% cvx_end
% 
% x_cvx = x;
% lambda_cvx = lambda;
% 
% 
% 
% norm(lambda_cvx,'inf')
% 
% mul_SASC = 1;
% eta_SASC = 2;
% 
% mul_SGD = 1.5;
% eta_SGD = 2;
% 
% mul_SGDM = 1.5;
% eta_SGDM = 2;
% alpha_SGDM = 0.9;
% 
% mul_SVRG = 1;
% eta_SVRG = 2;
% 
% mul_SVRG_screening = 1;
% eta_SVRG_screening = 2;
% 
% SVRG_ratio = 3;
% 
% Itrmax = 1e7 + 1;
% 
% sample = unidrnd(m,Itrmax,1); 
% 
% L = w + norm(Phi,2) ^ 2 / l;
% 
% mu = w;
% 
% err_SASC_hist = zeros(Itrmax,1);
% err_nested_SGD_hist = zeros(Itrmax,1);
% err_nested_SGDM_hist = zeros(Itrmax,1);
% err_static_SGDM_hist = zeros(Itrmax,1);
% err_nested_SVRG_hist = zeros(Itrmax,1);
% err_nested_SVRG_screening_hist = zeros(Itrmax,1);
% 
% 
% % SASC-SGD
% 
% x = x0;
% mul = mul_SASC;
% eta = eta_SASC;
% itr = 0;
% 
% inner_num = 0;
% 
% s = 1 / 6 / L; 
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
% 
%     beta = 4 * s;
% 
%     M = ceil(mul * eta / mu / s );
% 
%     for i = 1 : M
%         itr = itr + 1;
%         j = sample(itr);
%         
%         g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + (A(j,:))'*max(A(j,:)*x-b(j),0) / beta;
% 
%         x = x - s * g;
% 
%         err_SASC_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax
%         break
%     end
% 
%     s = s / eta;
% 
% end
% 
% inner_num_SASC = inner_num;
% 
% % Nested-SGD
% 
% x = x0;
% mul = mul_SGD;
% eta = eta_SGD;
% itr = 0;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
%     
%     s = 2 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));
% 
%     for i = 1 : M
%         itr = itr + 1;
%         j = sample(itr);
%         
%         g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
% 
%         x = x - s * g;
% 
%         err_nested_SGD_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax
%         break
%     end
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_SGD = inner_num;
% 
% % Nested-SGDM
% 
% x = x0;
% 
% v = zeros(n,1);
% 
% mul = mul_SGDM;
% eta = eta_SGDM;
% alpha = alpha_SGDM;
% itr = 0;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
%     
%     s = 2 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + m * xi / 4 / mu / delta ) ) * sqrt( m * ( L / mu + m * xi / 4 / mu / delta ) ) );
%     
%     v = zeros(n,1);
% 
%     for i = 1 : M
%         itr = itr + 1;
%         j = sample(itr);
%         
%         g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
%         
%         v = alpha * v - s * g;
% 
%         x = x + v;
% 
%         err_nested_SGDM_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax
%         break
%     end
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_SGDM = inner_num;
% 
% delta_nested_SGDM = delta;
% 
% % Static-SGDM
% 
% x = x0;
% 
% v = zeros(n,1);
% 
% delta = sqrt(delta_nested_SGDM *delta0);
% 
% s = 2 / ( L + mu + m * xi / 4 / delta );
% 
% for itr = 1 : Itrmax
%     j = sample(itr);
% 
%     g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
% 
%     v = alpha * v - s * g;
% 
%     x = x + v;
% 
%     err_static_SGDM_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% end
% 
% 
% SPDHG

err_SPDHG_hist = zeros(Itrmax,1);

% tau = 1 /  LA;
% 
% sigma = 1 / LA ;

gamma = 0.99;

sigma = ones(m,2);

for i = 1 : m
    sigma(i,2) = norm(Phi(i,:)) / sqrt(l);
end

tau = gamma / m / max(max(sigma));

sigma = gamma./ sigma;


x = x0;

u = zeros(m,1);

v = zeros(l,1);
f = A' * u;

f1 = f;

h = Phi' * v;

h1 = h;

for i = 1 : Itrmax
    x = ( x - tau * f1 - tau * h1) / ( w * tau + 1 );


    j = sample(i);
    utemp = u(j) + sigma(j,1) * A(j,:) * x;
    utemp2 = max(0,utemp-sigma(j,1)*b(j));
    
    vtemp = (v(j)+sigma(j,2)/sqrt(l)*(Phi(j,:)*x-y(j))) / (1+sigma(j,2));

    theta = 1/sqrt(1+2*mu*tau);
    tau = theta * tau;
    sigma = sigma / theta;


    f = f + (A(j,:))' * (utemp2-u(j));

    f1 = f + (A(j,:))' * (utemp2-u(j)) * m * theta;
    
    h = h + (Phi(j,:))'/sqrt(l) * (vtemp-v(j));

    h1 = h + (Phi(j,:))'/sqrt(l) * (vtemp-v(j)) * l * theta;

    u(j) = utemp2;
    
    v(j) = vtemp;

    err_SPDHG_hist(i) = norm(x-x_cvx,2)/norm(x_cvx,2); 
end
% 
% Itrmax1 = Itrmax;
% Itrmax2 = Itrmax;
% 
% sample2 = unidrnd(m,Itrmax2,1);
% 
% % Nested-SVRG
% 
% x = x0;
% mul = mul_SVRG;
% eta = eta_SVRG;
% itr = 0;
% N = SVRG_ratio;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
%     
%     s = 3 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));
% 
%     M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m; 
% 
%     i = 0;
% 
%     while i < M
%         if mod(i,(N+1)*m) == 0
%             z = x;
%             moment1 = zeros(n,1);
%             for j = 1 : m
%                 moment1 = moment1 + ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) / l + xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
%             end
%             i = i + m;
%         end
% 
%         itr = itr + 1;
%         i = i + 1;
% 
%         j = sample2(itr);
%         
%         gx = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
%         gz = ( Phi(j,:) )' * ( Phi(j,:) * z - y(j) ) + m * xi * ( A(j,:) )' * expd((A(j,:)*z-b(j))/delta);
%         
%         gdiff = gx - gz;
% 
%         x = ( x - s * gdiff - s * moment1 ) / ( 1 + s * w );
% 
%         err_nested_SVRG_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax2
%         break
%     end
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_SVRG = inner_num;
% 
% % Nested-SVRG-Screening
% 
% x = x0;
% mul = mul_SVRG_screening;
% eta = eta_SVRG_screening;
% itr = 0;
% N = SVRG_ratio;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% mtemp = m;
% Itemp = 1 : m;
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
%     
%     s = 3 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + mtemp * xi / 4 / mu / delta));
% 
%     M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m; 
% 
%     sample_temp = unidrnd(mtemp,M,1);
% 
%     i = 0;
% 
%     obj_i = N * m;
% 
%     cons_i = N * mtemp;
%     cons_j = 0;
% 
%     while i < 2 * M
%         if obj_i >= N * m
%             z = x;
%             moment1 = zeros(n,1);
%             for j = 1 : m
%                 moment1 = moment1 + ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) / l;
%             end
%             i = i + m;
% 
%             obj_i = 0;
%         end
%         if cons_i >= N * mtemp
%             v = x;
%             moment2 = zeros(n,1);
%             for k = 1 : mtemp
%                 j = Itemp(k);
%                 moment2 = moment2 + xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
%             end
%             i = i + mtemp;
% 
%             cons_i = 0;
%         end
%         itr = itr + 1;
%         obj_i = obj_i + 1;
%         cons_i = cons_i + 1;
%         cons_j = cons_j + 1;
%         i = i + 2;
% 
%         j = sample2(itr);
%         
%         gx1 = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) ;
%         gz1 = ( Phi(j,:) )' * ( Phi(j,:) * z - y(j) ) ;
%         
%         k = sample_temp(cons_j);
% 
%         gx2 = mtemp * xi * ( A(k,:) )' * expd((A(k,:)*x-b(k))/delta);
%         gv2 = mtemp * xi * ( A(k,:) )' * expd((A(k,:)*v-b(k))/delta);
% 
%         gdiff1 = gx1 - gz1;
%         gdiff2 = gx2 - gv2;
% 
%         x = ( x - s * gdiff1 - s * moment1 - s * gdiff2 - s * moment2 ) / ( 1 + s * w );
% 
%         err_nested_SVRG_screening_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax2
%         break
%     end
% 
%     Itemp_new = [];
% 
%     for i = 1 : mtemp
%         k = Itemp(i);
% 
%         val = A(k,:)*x-b(k);
% 
%         if val >= -2 * sqrt(mtemp) * delta
%             Itemp_new = [Itemp_new,k];
%         end
%     end
% 
%     Itemp = Itemp_new;
%     mtemp = length(Itemp);
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_SVRG_screening = inner_num;

% % semilogy(1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
% %     1:100:Itrmax1,err_nested_SGDM_hist(1:100:Itrmax1),1:100:Itrmax1,err_static_SGDM_hist(1:100:Itrmax1),...
% %     1:100:Itrmax1,err_SPDHG_hist(1:100:Itrmax1),(1+1/SVRG_ratio)*(1:100:Itrmax2),err_nested_SVRG_hist(1:100:Itrmax2));
% % legend('Nested-SGD','Nested-SGDM','Static-SGDM','SPDHG','Nested-SVRG');
% 
% semilogy(1:100:Itrmax,err_nested_SGD_hist(1:100:Itrmax),...
%     1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax),...
%     1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax),(1+1/SVRG_ratio)*(1:100:Itrmax),err_nested_SVRG_hist(1:100:Itrmax),...
%     (1+1/SVRG_ratio)*(1:100:Itrmax),err_nested_SVRG_screening_hist(1:100:Itrmax));
% legend('Nested-SGD','Nested-SGDM','Static-SGDM','SPDHG','Nested-SVRG','Nested-SVRG-Screening');

% semilogy(1:100:Itrmax,err_SASC_hist(1:100:Itrmax),1:100:Itrmax,err_nested_SGD_hist(1:100:Itrmax),...
%     1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax),...
%     1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax));
% legend('SASC-SGD','Nested-SGD','Nested-SGDM','Static-SGDM','SPDHG');


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


