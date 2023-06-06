% load("mushrooms_matrix.mat");
% load("mushrooms_vector.mat");
% 
% m = 8124;
% n = 111;
% 
% A = full(instance_matrix);
% 
% A = normalize(A);
% A = [A(:,1:77),A(:,79:112)];
% 
% b = -ones(m,1);
% 
% for i = 1 : m
%     nrm = norm(A(i,:),2);
%     A(i,:) = A(i,:) / nrm;
%     A(i,:) = - A(i,:) * label_vector(i);
%     b(i) = b(i) / nrm;
% end
% 
% rng(0)
% 
% 
% xi = 1;
% x0 = zeros(n,1);
% delta0 = 0.005;
% 
% cvx_begin
%     variable x(n,1);
%     dual variable lambda;
%     minimize 1 / 2 * (x'* x);
%     subject to
%         lambda: A * x <= b;
% cvx_end
% 
% x_cvx = x;
% lambda_cvx = lambda;
% 
% norm(lambda_cvx,'inf')
% 
% 
% Itrmax1 = 1.2e8 + 1;
% 
% Itrmax2 = 1e8 + 1;
% 
% sample = unidrnd(m,Itrmax1,1); 
% 
% L = 1;
% 
% mu = 1;
% 
% err_SASC_hist = zeros(Itrmax1,1);
% err_nested_SGD_hist = zeros(Itrmax1,1);
% err_nested_SVRG_hist = zeros(Itrmax1,1);
% err_nested_Screening_hist = zeros(Itrmax1,1);
% 
% 
% % SASC-SGD
% 
% x = x0;
% mul = 1;
% eta = 2;
% itr = 0;
% 
% inner_num = 0;
% 
% s = 1 / 2; 
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
%         g = (A(j,:))'*max(A(j,:)*x-b(j),0) / beta;
% 
%         x = ( x - s * g ) / ( 1 + s );
% 
%         err_SASC_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax1
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax1
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
% mul = 0.3;
% eta = 2;
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
%     s = 4 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));
% 
%     for i = 1 : M
%         itr = itr + 1;
%         j = sample(itr);
%         
%         g = m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);
% 
%         x = ( x - s * g ) / ( 1 + s );
% 
%         err_nested_SGD_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax1
%             break
%         end
% 
%     end
% 
%     if itr >= Itrmax1
%         break
%     end
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_SGD = inner_num;
% 
% % Nested-SVRG
%        
% sample2 = unidrnd(m,Itrmax2,1);
% x = x0;
% mul = 0.5;
% eta = 2;
% itr = 0;
% 
% N = 5;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
%     
%     s = 4 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));
% 
%     M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m; 
% 
%     i = 0;
% 
%     while i < M
%         if mod(i,(N+1)*m) == 0
%             y = x;
%             v = zeros(n,1);
%             for k = 1 : m
%                 v = v + xi * ( A(k,:) )' * expd((A(k,:)*x-b(k))/delta);
%             end
%             i = i + m;
%         end
% 
%         itr = itr + 1;
%         i = i + 1;
% 
%         j = sample2(itr);
%         
%         g = m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta) - m * xi * ( A(j,:) )' * expd((A(j,:)*y-b(j))/delta);
% 
%         x = ( x - s * g - s * v) / ( 1 + s );
% 
%         err_nested_SVRG_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%     end
% 
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
% % Nested-Screening
% 
% x = x0;
% mul = 0.5;
% eta = 2;
% itr = 0;
% 
% N = 5;
% 
% Itemp = 1:m;
% 
% mtemp = m;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     inner_num = inner_num + 1;
%     
%     s = 4 / ( L + mu + mtemp * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + mtemp * xi / 4 / mu / delta));
% 
%     M = ceil( M / ( N + 1 ) / mtemp ) * ( N + 1 ) * mtemp;
% 
%     sample_temp = unidrnd(mtemp,M,1);
%     
%     i = 0;
% 
%     while i < M
%         if mod(i,(N+1)*mtemp) == 0
%             y = x;
%             v = zeros(n,1);
%             for k = 1 : mtemp
%                 ki = Itemp(k);
%                 v = v + xi * ( A(ki,:) )' * expd((A(ki,:)*x-b(ki))/delta);
%             end
%             i = i + mtemp;
%         end
% 
%         itr = itr + 1;
%         i = i + 1;
% 
%         j = Itemp(sample_temp(i));
%         
%         g = mtemp * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta) - mtemp * xi * ( A(j,:) )' * expd((A(j,:)*y-b(j))/delta);
% 
%         x = ( x - s * g - s * v) / ( 1 + s );
% 
%         err_nested_Screening_hist(itr) = norm(x-x_cvx,2)/norm(x_cvx,2);
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%     end
% 
% 
%     if itr >= Itrmax2
%         break
%     end
%     
%     Itemp_new = [];
% 
%     for i = 1 : mtemp
%         ii = Itemp(i);
%         if A(ii,:) * x - b(ii) >= -2 * sqrt(mtemp) * delta
%             Itemp_new = [Itemp_new,ii];
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
% inner_num_Screening = inner_num;

% semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
%     1.2*(1:100:Itrmax2),err_nested_SVRG_hist(1:100:Itrmax2));
% legend('SASC-SGD','Nested-SGD','Nested-SVRG');
% xlabel('Number of Stochastic Gradient Evaluations');
% ylabel('Relative Error');

% semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
%     1.2*(1:100:Itrmax2),err_nested_SVRG_hist(1:100:Itrmax2),1.2*(1:100:Itrmax2),err_nested_Screening_hist(1:100:Itrmax2));
% legend('SASC-SGD','Nested-SGD','Nested-SVRG','Nested-SVRG-Screening');
% xlabel('Number of Stochastic Gradient Evaluations');
% ylabel('Relative Error');

log_Err = zeros(m,1);

log_Err_Itemp = zeros(mtemp,1);

for i = 1 : m
    Err_temp = b(i) - A(i,:) * x_cvx;
    if Err_temp <= 1e-10
        log_Err(i) = -10;
    else
        log_Err(i) = log10(Err_temp);
    end

end

for ij = 1 : mtemp
    i = Itemp(ij);
    Err_temp = b(i) - A(i,:) * x_cvx;
    if Err_temp <= 1e-10
        log_Err_Itemp(ij) = -10;
    else
        log_Err_Itemp(ij) = log10(Err_temp);
    end

end
edges = -10: 1: 2;
h1 = histogram(log_Err,edges);
hold on
h2 = histogram(log_Err_Itemp,edges);
legend('Original','Remaining')
xlabel('lg(Slackness)');
ylabel('Number of Constraints')
