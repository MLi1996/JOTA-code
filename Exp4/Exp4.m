% rng(3)
% 
% n = 500;
% l = 100;
% m = 100;
% 
% w_int = randn(n,1);
% 
% w_int = w_int / norm(w_int);
% 
% X10 = randn(l,n);
% Y1 = zeros(l,1);
% 
% for i = 1:l
%     X10(i,:) = X10(i,:) / norm(X10(i,:));
%     Y1(i) = ( X10(i,:) * w_int > 0 );
% end
% 
% X2 = randn(m,n);
% Y2 = zeros(m,1);
% 
% for i = 1:m
%     X2(i,:) = X2(i,:) / norm(X2(i,:));
%     Y2(i) = ( X2(i,:) * w_int > 0 );
% end
% 
% X11 = randn(l,n);
% 
% 
% for i = 1:l
%     X11(i,:) = X11(i,:) / norm(X11(i,:));
% end
% 
% X1 = X10 + 0.1 * X11;
% 
% 
% mu = 0.05;
% 
% 
% A2 = zeros( m , n );
% 
% for i = 1 : m
%     A2(i,:) = ( 2 * Y2(i) - 1 ) * X2(i,:);
% end 
% 
% 
% cvx_begin
%     variable w(n);
%     dual variable lambda;
%     minimize (- Y1' * X1 * w +sum(log_sum_exp([zeros(1,l); w'*X1'])) + mu / 2 * sum_square(w));
%     subject to
%         lambda: A2 * w >= 1;
% cvx_end
% 
% w_cvx = w;
% lambda_cvx = lambda;
% 
% 
% norm(lambda_cvx,'inf')
% 
% mul_SASC = 1;
% eta_SASC = 2;
% 
% mul_SGD = 1;
% eta_SGD = 2;
% 
% SVRG_ratio = 3;
% 
% mul_SVRG = 0.5;
% eta_SVRG = 2;
% 
% mul_SVRG_screening = 0.5;
% eta_SVRG_screening = 2;
% 
% mul_SGDM = 1;
% eta_SGDM = 2;
% alpha_SGDM = 0.5;
% 
% %%%%%%
% 
% 
% Itrmax = 1e6 + 1;
% 
% sample = unidrnd(m,Itrmax,1); 
% % sample2 = unidrnd(m,Itrmax,1);
% 
% num_C = 2;
% 
% S = kron((eye(num_C)-ones(num_C,num_C)/num_C), X1*X1') / 2;
% 
% L = max(eig(S)) + mu;
% 
% mu = mu;
% 
% err_SASC_hist = zeros(Itrmax,1);
% err_nested_SGD_hist = zeros(Itrmax,1);
% err_nested_SGDM_hist = zeros(Itrmax,1);
% err_static_SGDM_hist = zeros(Itrmax,1);
% err_nested_SVRG_hist = zeros(Itrmax,1);
% err_nested_SVRG_Screening_hist = zeros(Itrmax,1);
% 
% xi = 1;
% delta0 = 1;
% 
% 
% w0 = zeros(n,1);
% 
% % SASC-SGD
% 
% w = w0;
% mul = mul_SASC;
% eta = eta_SASC;
% itr = 0;
% 
% inner_num = 0;
% 
% s = 1 / 2 / L; 
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
%         g = -l * (Y1(j)*X1(j,:))';
% 
%         g = g + l * expd(X1(j,:)*w) * (X1(j,:))'; 
%         
%         val = max( -A2(j,:) * w + 1, 0 ) / beta;
%         g = g - val * (A2(j,:))';
% 
%         w = ( w - s * g ) / ( 1 + s * mu ) ;
% 
%         err_SASC_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
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
% w = w0;
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
%     s = 1 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));
% 
%     for i = 1 : M
%         
%         itr = itr + 1;
%         j = sample(itr);
%         
%         g = -l * (Y1(j)*X1(j,:))';
% 
%         g = g + l * expd(X1(j,:)*w) * (X1(j,:))'; 
%         
%         val = expd( (-A2(j,:) * w + 1) / delta ) * m * xi;
% 
%         g = g - val * (A2(j,:))';
% 
%         w = ( w - s * g ) / ( 1 + s * mu );
% 
%         err_nested_SGD_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
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
% 
% Itrmax1 = Itrmax;
% Itrmax2 = Itrmax;
% 
% % Nested-SVRG
%        
% sample2 = unidrnd(m,Itrmax2,1);
% w = w0;
% mul = mul_SVRG;
% eta = eta_SVRG;
% N = SVRG_ratio;
% 
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
%     M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m; 
% 
%     i = 0;
% 
%     while i < M
%         if mod(i,(N+1)*m) == 0
%             y = w;
%             v1 = zeros(n,1);
%             for j = 1 : m
% 
%                 g = - (Y1(j)*X1(j,:))';
% 
%                 g = g + expd(X1(j,:)*w) * (X1(j,:))';
% 
%                 val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
% 
%                 g = g - val * (A2(j,:))';
% 
%                 v1 = v1 + g;
% 
%             end
% 
%             i = i + m;
% 
%         end
% 
%         itr = itr + 1;
%         i = i + 1;
% 
%         j = sample2(itr);
%         
%         gx = - l * (Y1(j)*X1(j,:))';
% 
%         gx = gx + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%         val = expd( (-A2(j,:) * w + 1) / delta ) * m * xi;
% 
%         gx = gx - val * (A2(j,:))';
% 
%         gy = - l * (Y1(j)*X1(j,:))';
% 
%         gy = gy + l * expd(X1(j,:)*y) * (X1(j,:))';
% 
%         val = expd( (-A2(j,:) * y + 1) / delta ) * m * xi;
% 
%         gy = gy - val * (A2(j,:))';
% 
%         gdiff = gx - gy;
% 
%         w = ( w - s * gdiff - s * v1) / ( 1 + s * mu );
% 
%         err_nested_SVRG_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
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
% 
% 
% % Nested-SVRG_Screening
% 
% w = w0;
% mul = mul_SVRG_screening;
% eta = eta_SVRG_screening;
% N = SVRG_ratio;
% 
% mtemp = m;
% 
% Itemp = 1:mtemp;
% 
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
%     M = ceil(mul * log(2*eta-1) * ( L / mu + mtemp * xi / 4 / mu / delta));
% 
%     M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m;
% 
%     sample_temp = unidrnd(mtemp,M,1);
% 
%     i = 0;
% 
%     obj_i = N * l;
% 
%     cons_i = N * mtemp;
% 
%     cons_j = 0;
% 
%     while i < 2 * M
%         if obj_i >= N * l
%             y = w;
%             v1 = zeros(n,1);
%             for j = 1 : m
%                 g = - (Y1(j)*X1(j,:))';
%                 g = g + expd(X1(j,:)*w) * (X1(j,:))';
%                 v1 = v1 + g;
%             end
% 
%             i = i + m;
%             obj_i = 0;
%         end
%             
%         if cons_i >= N * mtemp
%             z = w;
%             v2 = zeros(n,1);
%             for ktemp = 1 : mtemp
%                 j = Itemp(ktemp);
%                 val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
%                 g = - val * (A2(j,:))';
%                 v2 = v2 + g;
%             end
%             i = i + mtemp;
%             cons_i = 0;
%         end
% 
% 
%         itr = itr + 1;
%         i = i + 2;
%         obj_i = obj_i + 1;
%         cons_i = cons_i + 1;
%         cons_j = cons_j + 1;
% 
%         j_cons = Itemp(sample_temp(cons_j));
%         j = sample2(itr);
% 
%         gx_obj = - l * (Y1(j)*X1(j,:))';
% 
%         gx_obj = gx_obj + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%         gy_obj = - l * (Y1(j)*X1(j,:))';
% 
%         gy_obj = gy_obj + l * expd(X1(j,:)*y) * (X1(j,:))';
% 
%         gdiff_obj = gx_obj - gy_obj;
% 
%         val = expd( (-A2(j_cons,:) * w + 1) / delta ) * xi * mtemp ;
% 
%         gx_cons = - val * (A2(j_cons,:))';
% 
%         val = expd( (-A2(j_cons,:) * z + 1) / delta ) * xi * mtemp ;
% 
%         gz_cons = - val * (A2(j_cons,:))';
% 
%         gdiff_cons = gx_cons - gz_cons;
% 
%         w = ( w - s * gdiff_obj - s * v1 - s * gdiff_cons - s * v2 ) / ( 1 + s * mu );
% 
%         err_nested_SVRG_Screening_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%     end
% 
% 
% 
%     if itr >= Itrmax2
%         break
%     end
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
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_Screening = inner_num;
% 
% 
% 
% % Nested-SGDM
% 
% w = w0;
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
%     s = 1 / ( L + mu + m * xi / 4 / delta );
% 
%     M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + m * xi / 4 / mu / delta ) ) * sqrt( m * ( L / mu + m * xi / 4 / mu / delta ) ) );
%     
%     v = zeros(n,1);
% 
%     for i = 1 : M
%         itr = itr + 1;
%         j = sample(itr);
%         
%         g = -l * (Y1(j)*X1(j,:))';
% 
%         g = g + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%         val = expd( (-A2(j,:) * w + 1) / delta ) * m * xi;
% 
%         g = g - val * (A2(j,:))' + mu * w;
% 
%         v = alpha * v - s * g;
% 
%         w = w + v;
% 
%         err_nested_SGDM_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
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
% w = w0;
% 
% v = zeros(n,1);
% 
% delta = sqrt(delta_nested_SGDM *delta0);
% 
% s = 1 / ( L + mu + m * xi / 4 / delta );
% 
% for itr = 1 : Itrmax
%     j = sample(itr);
% 
%     g = -l * (Y1(j)*X1(j,:))';
% 
%     g = g + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%     val = expd( (-A2(j,:) * w + 1) / delta ) * m * xi;
% 
%     g = g - val * (A2(j,:))' + mu * w;
% 
%     v = alpha * v - s * g;
% 
%     w = w + v;
% 
%     err_static_SGDM_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
% end
% 
% 
% % SPDHG
% 
% err_SPDHG_hist = zeros(Itrmax,1);
% 
% % tau = 1 /  LA;
% % 
% % sigma = 1 / LA ;
% 
% 
% Sigma = ones(m,2);
% 
% for i = 1 : m
%     Sigma(i,2) = norm(X1(i,:),2);
% end
% 
% gamma = 0.99;
% 
% tau = gamma / m / max(max(Sigma)) / 2;
% 
% sigma = gamma./ Sigma / 2 ;
% 
% 
% w = w0;
% 
% psi = zeros(n,1);
% 
% for i = 1 : l
%     psi = psi + Y1(i) * (X1(i,:))'; 
% end
% 
% u = zeros(m,1);
% 
% v = zeros(l,1);
% 
% f = A2' * u;
% 
% f1 = f;
% 
% h = X1' * v;
% 
% h1 = h;
% 
% for i = 1 : Itrmax
%     w = ( w - tau * f1 - tau * h1) / ( mu * tau + 1 ) +  psi * tau / ( 1 + tau * mu ) ;
% 
% 
%     j = sample(i);
%     utemp = min(0, u(j) + sigma(j,1) * (A2(j,:) * w - 1) );
% 
%     
%     vtemp0 = v(j) + sigma(j,2) * X1(j,:) * w;
% 
%     func = @(val) sigma(j,2) * ( log(val) - log(1-val) ) + ( val - vtemp0 );
% 
%     options.Display = 'off';
%     options.MaxIterations = 5;
% 
%     init = v(j);
% 
%     if init < 1e-10 || init > 1-1e-10
%         init = 0.5;
%     end
% 
% 
%     vtemp = real(fsolve(func,init,options));
%     
%  
%     theta = 1/sqrt(1+2*mu*tau);
%     tau = theta * tau;
%     sigma = sigma / theta;
% 
% 
%     f = f + (A2(j,:))' * (utemp-u(j));
% 
%     f1 = f + (A2(j,:))' * (utemp-u(j)) * m * theta;
%     
%     h = h + (X1(j,:))' * (vtemp-v(j));
% 
%     h1 = h + (X1(j,:))' * (vtemp-v(j)) * l * theta;
% 
%     u(j) = utemp;
%     
%     v(j) = vtemp;
% 
%     err_SPDHG_hist(i) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
% end


% 
% % Itrmax1 = Itrmax;
% % 
% % semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
% %      1.33333*(1:100:Itrmax2),err_nested_SVRG_hist(1:100:Itrmax2),...
% %      1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax));
% % legend('SASC-SGD','Nested-SGD','Nested-SVRG','Nested-SGDM','Static-SGDM');
% % xlabel('Number of Stochastic Evaluations');
% % ylabel('Relative Error');
% % 
% % 
% semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
%      1.33333*(1:100:Itrmax2),err_nested_SVRG_hist(1:100:Itrmax2),1.33333*(1:100:Itrmax2),err_nested_SVRG_Screening_hist(1:100:Itrmax2),...
%      1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax));
% legend('SASC-SGD','Nested-SGD','Nested-SVRG','Nested-SVRG-Screening','Nested-SGDM','Static-SGDM');
% xlabel('Number of Stochastic Evaluations');
% ylabel('Relative Error');


% semilogy(1:100:Itrmax1,err_SASC_hist(1:100:Itrmax1),1:100:Itrmax1,err_nested_SGD_hist(1:100:Itrmax1),...
%      1.33333*(1:100:Itrmax2),err_nested_SVRG_hist(1:100:Itrmax2),1.33333*(1:100:Itrmax2),err_nested_SVRG_Screening_hist(1:100:Itrmax2),...
%      1:100:Itrmax,err_nested_SGDM_hist(1:100:Itrmax),1:100:Itrmax,err_static_SGDM_hist(1:100:Itrmax),...
%      1:100:Itrmax,err_SPDHG_hist(1:100:Itrmax));
% legend('SASC-SGD','Nested-SGD','Nested-SVRG','Nested-SVRG-Screening','Nested-SGDM','Static-SGDM','SPDHG');
% xlabel('Number of Stochastic Evaluations');
% ylabel('Relative Error');



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

% save('Exp4_SVRG_Screening.mat','err_nested_SVRG_Screening_hist');
% save('Exp4_SPDHG.mat','err_SPDHG_hist');

