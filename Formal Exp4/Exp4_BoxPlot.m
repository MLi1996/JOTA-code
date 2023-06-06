% rng(3)
% 
% n = 500;
% l = 100;
% m = 100;
% 
% 
% Exp_Time = 10;
% 
% 
% Itrmax = 1e6 + 1;
% 
% 
% err_nested_SVRG_hist = zeros(Itrmax,1);
% err_nested_SVRG_Screening_hist = zeros(Itrmax,1);
% err_SPDHG_hist = zeros(Itrmax,1);
% err_nested_Catalyst_hist = zeros(Itrmax,1);
% err_nested_Catalyst_Screening_hist = zeros(Itrmax,1);
% 
% 
% Err_SVRG = zeros(Exp_Time,1);
% Err_SVRG_Screening = zeros(Exp_Time,1);
% Err_SPDHG = zeros(Exp_Time,1);
% Err_Catalyst = zeros(Exp_Time,1);
% Err_Catalyst_Screening = zeros(Exp_Time,1);
% 
% 
% for exp_time = 1 : Exp_Time
%     w_int = randn(n,1);
% 
%     w_int = w_int / norm(w_int);
% 
%     X10 = randn(l,n);
%     Y1 = zeros(l,1);
% 
%     for i = 1:l
%         X10(i,:) = X10(i,:) / norm(X10(i,:));
%         Y1(i) = ( X10(i,:) * w_int > 0 );
%     end
% 
%     X2 = randn(m,n);
%     Y2 = zeros(m,1);
% 
%     for i = 1:m
%         X2(i,:) = X2(i,:) / norm(X2(i,:));
%         Y2(i) = ( X2(i,:) * w_int > 0 );
%     end
% 
%     X11 = randn(l,n);
% 
% 
%     for i = 1:l
%         X11(i,:) = X11(i,:) / norm(X11(i,:));
%     end
% 
%     X1 = X10 + 0.1 * X11;
% 
% 
%     mu = 0.05;
% 
% 
%     A2 = zeros( m , n );
% 
%     for i = 1 : m
%         A2(i,:) = ( 2 * Y2(i) - 1 ) * X2(i,:);
%     end
% 
% 
%     cvx_begin
%     variable w(n);
%     dual variable lambda;
%     minimize (- Y1' * X1 * w +sum(log_sum_exp([zeros(1,l); w'*X1'])) + mu / 2 * sum_square(w));
%     subject to
%     lambda: A2 * w >= 1;
%     cvx_end
% 
%     w_cvx = w;
%     lambda_cvx = lambda;
% 
% 
%     norm(lambda_cvx,'inf')
% 
%     mul_SASC = 1;
%     eta_SASC = 2;
% 
%     mul_SGD = 1;
%     eta_SGD = 2;
% 
%     SVRG_ratio = 3;
%     Catalyst_Budget = 1;
% 
%     mul_SVRG = 0.5;
%     eta_SVRG = 2;
% 
%     mul_SVRG_screening = 0.5;
%     eta_SVRG_screening = 2;
% 
%     mul_Catalyst = 1.5;
%     eta_Catalyst = 2;
% 
%     mul_Catalyst_Screening = 1.5;
%     eta_Catalyst_Screening = 2;
% 
%     mul_SGDM = 1;
%     eta_SGDM = 2;
%     alpha_SGDM = 0.5;
% 
%     %%%%%%
% 
% 
%     num_C = 2;
% 
%     S = kron((eye(num_C)-ones(num_C,num_C)/num_C), X1*X1') / 2;
% 
%     L = max(eig(S)) + mu;
% 
%     mu = mu;
% 
% 
%     xi = 1;
%     delta0 = 1;
% 
% 
%     w0 = zeros(n,1);
% 
% 
%     Itrmax1 = Itrmax;
%     Itrmax2 = Itrmax;
%     sample = unidrnd(m,Itrmax1,1);
% 
%     % Nested-SVRG
% 
%     sample2 = unidrnd(m,Itrmax2,1);
%     w = w0;
%     mul = mul_SVRG;
%     eta = eta_SVRG;
%     N = SVRG_ratio;
% 
%     itr = 0;
% 
%     inner_num = 0;
% 
%     delta = delta0;
% 
%     while inner_num >= 0
% 
%         inner_num = inner_num + 1;
% 
%         s = 2 / ( L + mu + m * xi / 4 / delta );
% 
%         M = ceil(mul * log(2*eta-1) * ( L / mu + m * xi / 4 / mu / delta));
% 
%         M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m;
% 
%         i = 0;
% 
%         while i < M
%             if mod(i,(N+1)*m) == 0
%                 y = w;
%                 v1 = zeros(n,1);
%                 for j = 1 : m
% 
%                     g = - (Y1(j)*X1(j,:))';
% 
%                     g = g + expd(X1(j,:)*w) * (X1(j,:))';
% 
%                     val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
% 
%                     g = g - val * (A2(j,:))';
% 
%                     v1 = v1 + g;
% 
%                 end
% 
%                 i = i + m;
% 
%             end
% 
%             itr = itr + 1;
%             i = i + 1;
% 
%             j = sample2(itr);
% 
%             gx = - l * (Y1(j)*X1(j,:))';
% 
%             gx = gx + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%             val = expd( (-A2(j,:) * w + 1) / delta ) * m * xi;
% 
%             gx = gx - val * (A2(j,:))';
% 
%             gy = - l * (Y1(j)*X1(j,:))';
% 
%             gy = gy + l * expd(X1(j,:)*y) * (X1(j,:))';
% 
%             val = expd( (-A2(j,:) * y + 1) / delta ) * m * xi;
% 
%             gy = gy - val * (A2(j,:))';
% 
%             gdiff = gx - gy;
% 
%             w = ( w - s * gdiff - s * v1) / ( 1 + s * mu );
% 
%             err_nested_SVRG_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
%             if itr >= Itrmax2
%                 break
%             end
% 
%         end
% 
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%         delta = delta / eta;
% 
%     end
% 
%     inner_num_SVRG = inner_num;
% 
% 
% 
%     % Nested-SVRG_Screening
% 
%     w = w0;
%     mul = mul_SVRG_screening;
%     eta = eta_SVRG_screening;
%     N = SVRG_ratio;
% 
%     mtemp = m;
% 
%     Itemp = 1:mtemp;
% 
%     itr = 0;
% 
%     inner_num = 0;
% 
%     delta = delta0;
% 
%     while inner_num >= 0
% 
%         inner_num = inner_num + 1;
% 
%         s = 2 / ( L + mu + m * xi / 4 / delta );
% 
%         M = ceil(mul * log(2*eta-1) * ( L / mu + mtemp * xi / 4 / mu / delta));
% 
%         M = ceil( M / ( N + 1 ) / m ) * ( N + 1 ) * m;
% 
%         sample_temp = unidrnd(mtemp,M,1);
% 
%         i = 0;
% 
%         obj_i = N * l;
% 
%         cons_i = N * mtemp;
% 
%         cons_j = 0;
% 
%         while i < 2 * M
%             if obj_i >= N * l
%                 y = w;
%                 v1 = zeros(n,1);
%                 for j = 1 : m
%                     g = - (Y1(j)*X1(j,:))';
%                     g = g + expd(X1(j,:)*w) * (X1(j,:))';
%                     v1 = v1 + g;
%                 end
% 
%                 i = i + m;
%                 obj_i = 0;
%             end
% 
%             if cons_i >= N * mtemp
%                 z = w;
%                 v2 = zeros(n,1);
%                 for ktemp = 1 : mtemp
%                     j = Itemp(ktemp);
%                     val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
%                     g = - val * (A2(j,:))';
%                     v2 = v2 + g;
%                 end
%                 i = i + mtemp;
%                 cons_i = 0;
%             end
% 
% 
%             itr = itr + 1;
%             i = i + 2;
%             obj_i = obj_i + 1;
%             cons_i = cons_i + 1;
%             cons_j = cons_j + 1;
% 
%             j_cons = Itemp(sample_temp(cons_j));
%             j = sample2(itr);
% 
%             gx_obj = - l * (Y1(j)*X1(j,:))';
% 
%             gx_obj = gx_obj + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%             gy_obj = - l * (Y1(j)*X1(j,:))';
% 
%             gy_obj = gy_obj + l * expd(X1(j,:)*y) * (X1(j,:))';
% 
%             gdiff_obj = gx_obj - gy_obj;
% 
%             val = expd( (-A2(j_cons,:) * w + 1) / delta ) * xi * mtemp ;
% 
%             gx_cons = - val * (A2(j_cons,:))';
% 
%             val = expd( (-A2(j_cons,:) * z + 1) / delta ) * xi * mtemp ;
% 
%             gz_cons = - val * (A2(j_cons,:))';
% 
%             gdiff_cons = gx_cons - gz_cons;
% 
%             w = ( w - s * gdiff_obj - s * v1 - s * gdiff_cons - s * v2 ) / ( 1 + s * mu );
% 
%             err_nested_SVRG_Screening_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
%             if itr >= Itrmax2
%                 break
%             end
% 
%         end
% 
% 
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%         Itemp_new = [];
% 
%         for i = 1 : mtemp
% 
%             j = Itemp(i);
% 
%             val = -A2(j,:) * w + 1;
% 
%             if val >= -2 * sqrt(mtemp) * delta
%                 Itemp_new = [Itemp_new,j];
%             end
%         end
% 
%         Itemp = Itemp_new;
%         mtemp = length(Itemp);
% 
%         delta = delta / eta;
% 
%     end
% 
%     inner_num_Screening = inner_num;
% 
%     % SPDHG
% 
%     Sigma = ones(m,2);
% 
%     for i = 1 : m
%         Sigma(i,2) = norm(X1(i,:),2);
%     end
% 
%     gamma = 0.99;
% 
%     tau = gamma / m / max(max(Sigma)) / 2;
% 
%     sigma = gamma./ Sigma / 2 ;
% 
% 
%     w = w0;
% 
%     psi = zeros(n,1);
% 
%     for i = 1 : l
%         psi = psi + Y1(i) * (X1(i,:))';
%     end
% 
%     u = zeros(m,1);
% 
%     v = zeros(l,1);
% 
%     f = A2' * u;
% 
%     f1 = f;
% 
%     h = X1' * v;
% 
%     h1 = h;
% 
%     for i = 1 : Itrmax1
%         w = ( w - tau * f1 - tau * h1) / ( mu * tau + 1 ) +  psi * tau / ( 1 + tau * mu ) ;
% 
% 
%         j = sample(i);
%         utemp = min(0, u(j) + sigma(j,1) * (A2(j,:) * w - 1) );
% 
% 
%         vtemp0 = v(j) + sigma(j,2) * X1(j,:) * w;
% 
%         func = @(val) sigma(j,2) * ( log(val) - log(1-val) ) + ( val - vtemp0 );
% 
%         options.Display = 'off';
%         options.MaxIterations = 5;
% 
%         init = v(j);
% 
%         if init < 1e-10 || init > 1-1e-10
%             init = 0.5;
%         end
% 
% 
%         vtemp = real(fsolve(func,init,options));
% 
% 
%         theta = 1/sqrt(1+2*mu*tau);
%         tau = theta * tau;
%         sigma = sigma / theta;
% 
% 
%         f = f + (A2(j,:))' * (utemp-u(j));
% 
%         f1 = f + (A2(j,:))' * (utemp-u(j)) * m * theta;
% 
%         h = h + (X1(j,:))' * (vtemp-v(j));
% 
%         h1 = h + (X1(j,:))' * (vtemp-v(j)) * l * theta;
% 
%         u(j) = utemp;
% 
%         v(j) = vtemp;
% 
%         err_SPDHG_hist(i) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
%     end
% 
% 
% 
%     % Nested-Catalyst_SVRG
% 
%     ritr = 0;
% 
%     w0 = zeros(n,1);
% 
%     w = w0;
%     mul = mul_Catalyst;
%     eta = eta_Catalyst;
%     N = SVRG_ratio;
% 
%     T_B = Catalyst_Budget;
% 
%     mtemp = m;
% 
%     Itemp = 1:mtemp;
% 
%     itr = 0;
% 
%     inner_num = 0;
% 
%     delta = delta0;
% 
%     while inner_num >= 0
% 
%         wy = w;
% 
%         wold = w;
% 
%         inner_num = inner_num + 1;
% 
%         Ltemp = L + m * xi / 4 / delta;
% 
%         kappa = ( L - mu ) / m ;
% 
%         q = mu / ( mu + kappa );
% 
%         alpha = sqrt(q);
% 
%         beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );
% 
%         s = 3 / ( Ltemp + mu + kappa );
% 
%         M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + mtemp * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + mtemp^2 * xi / 4 / mu / delta  ) );
% 
%         M = ceil( M / ( N + 1 ) / mtemp / T_B) * ( N + 1 ) * mtemp * T_B;
% 
%         sample_temp = unidrnd(mtemp,M,1);
% 
%         i = 0;
%         t_b = 0;
% 
%         while i < M
% 
%             if mod(i,(N+1)*mtemp) == 0
%                 y = w;
%                 v1 = zeros(n,1);
%                 for j = 1 : m
% 
%                     g = - (Y1(j)*X1(j,:))';
%                     g = g + expd(X1(j,:)*w) * (X1(j,:))';
%                     v1 = v1 + g;
% 
%                 end
%                 ritr = ritr + m;
% 
%                 for jj = 1 : mtemp
% 
%                     j = Itemp(jj);
%                     val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
%                     v1 = v1 - val * (A2(j,:))';
% 
%                 end
% 
%                 ritr = ritr + mtemp;
%                 i = i + mtemp;
% 
%             end
% 
%             ritr = ritr + 2;
% 
%             itr = itr + 1;
%             i = i + 1;
% 
% 
%             j = sample2(itr);
%             j_cons = Itemp(sample_temp(i));
% 
%             gx = - l * (Y1(j)*X1(j,:))';
% 
%             gx = gx + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%             val = expd( (-A2(j_cons,:) * w + 1) / delta ) * mtemp * xi;
% 
%             gx = gx - val * (A2(j_cons,:))';
% 
%             gy = - l * (Y1(j)*X1(j,:))';
% 
%             gy = gy + l * expd(X1(j,:)*y) * (X1(j,:))';
% 
%             val = expd( (-A2(j_cons,:) * y + 1) / delta ) * mtemp * xi;
% 
%             gy = gy - val * (A2(j_cons,:))';
% 
%             gdiff = gx - gy;
% 
%             w = ( w - s * gdiff - s * v1) / ( 1 + s * mu + s * kappa ) + wy * s * kappa / ( 1 + s * mu + s * kappa );
% 
%             err_nested_Catalyst_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
%             if itr >= Itrmax2
%                 break
%             end
% 
%             if mod(i,(N+1)*mtemp) == 0
%                 t_b = t_b + 1;
%             end
% 
%             if t_b == T_B
%                 wy = w + beta * (w - wold);
%                 wold = w;
%                 t_b = 0;
%             end
%         end
% 
% 
% 
%         if itr >= Itrmax2
%             break;
%         end
%         %
%         %     Itemp_new = [];
%         %
%         %     for i = 1 : mtemp
%         %
%         %         j = Itemp(i);
%         %
%         %         val = -A2(j,:) * w + 1;
%         %
%         %         if val >= -2 * sqrt(mtemp) * delta
%         %             Itemp_new = [Itemp_new,j];
%         %         end
%         %     end
%         %
%         %     Itemp = Itemp_new;
%         %     mtemp = length(Itemp);
% 
%         delta = delta / eta;
% 
%     end
% 
%     ritr_Catalyst = ritr / 2;
% 
%     inner_num_catalyst = inner_num;
% 
% 
% 
% 
% 
% 
%     % Nested-Catalyst_SVRG_Screening
% 
%     ritr = 0;
% 
%     w0 = zeros(n,1);
% 
%     w = w0;
%     mul = mul_Catalyst_Screening;
%     eta = eta_Catalyst_Screening;
%     N = SVRG_ratio;
% 
%     T_B = Catalyst_Budget;
% 
%     mtemp = m;
% 
%     Itemp = 1:mtemp;
% 
%     itr = 0;
% 
%     inner_num = 0;
% 
%     delta = delta0;
% 
%     while inner_num >= 0
% 
%         wy = w;
% 
%         wold = w;
% 
%         inner_num = inner_num + 1;
% 
%         Ltemp = L + m * xi / 4 / delta;
% 
%         kappa = ( L - mu ) / m ;
% 
%         q = mu / ( mu + kappa );
% 
%         alpha = sqrt(q);
% 
%         beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );
% 
%         s = 3 / ( Ltemp + mu + kappa );
% 
%         M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + mtemp * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + mtemp^2 * xi / 4 / mu / delta  ) );
% 
%         M = ceil( M / ( N + 1 ) / mtemp / T_B) * ( N + 1 ) * mtemp * T_B;
% 
%         sample_temp = unidrnd(mtemp,M,1);
% 
%         i = 0;
%         t_b = 0;
% 
%         while i < M
% 
%             if mod(i,(N+1)*mtemp) == 0
%                 y = w;
%                 v1 = zeros(n,1);
%                 for j = 1 : m
% 
%                     g = - (Y1(j)*X1(j,:))';
%                     g = g + expd(X1(j,:)*w) * (X1(j,:))';
%                     v1 = v1 + g;
% 
%                 end
% 
%                 ritr = ritr + m;
% 
%                 for jj = 1 : mtemp
% 
%                     j = Itemp(jj);
%                     val = expd( (-A2(j,:) * w + 1) / delta ) * xi;
%                     v1 = v1 - val * (A2(j,:))';
% 
%                 end
% 
%                 ritr = ritr + mtemp;
% 
%                 i = i + mtemp;
% 
%             end
% 
%             ritr = ritr + 2;
% 
%             itr = itr + 1;
%             i = i + 1;
% 
% 
%             j = sample2(itr);
%             j_cons = Itemp(sample_temp(i));
% 
%             gx = - l * (Y1(j)*X1(j,:))';
% 
%             gx = gx + l * expd(X1(j,:)*w) * (X1(j,:))';
% 
%             val = expd( (-A2(j_cons,:) * w + 1) / delta ) * mtemp * xi;
% 
%             gx = gx - val * (A2(j_cons,:))';
% 
%             gy = - l * (Y1(j)*X1(j,:))';
% 
%             gy = gy + l * expd(X1(j,:)*y) * (X1(j,:))';
% 
%             val = expd( (-A2(j_cons,:) * y + 1) / delta ) * mtemp * xi;
% 
%             gy = gy - val * (A2(j_cons,:))';
% 
%             gdiff = gx - gy;
% 
%             w = ( w - s * gdiff - s * v1) / ( 1 + s * mu + s * kappa ) + wy * s * kappa / ( 1 + s * mu + s * kappa );
% 
%             err_nested_Catalyst_Screening_hist(itr) = norm(w-w_cvx,2)/norm(w_cvx,2);
% 
%             if itr >= Itrmax2
%                 break
%             end
% 
%             if mod(i,(N+1)*mtemp) == 0
%                 t_b = t_b + 1;
%             end
% 
%             if t_b == T_B
%                 wy = w + beta * (w - wold);
%                 wold = w;
%                 t_b = 0;
%             end
%         end
% 
% 
% 
%         if itr >= Itrmax2
%             break;
%         end
% 
%         Itemp_new = [];
% 
%         for i = 1 : mtemp
% 
%             j = Itemp(i);
% 
%             val = -A2(j,:) * w + 1;
% 
%             if val >= -2 * sqrt(mtemp) * delta
%                 Itemp_new = [Itemp_new,j];
%             end
%         end
% 
%         Itemp = Itemp_new;
%         mtemp = length(Itemp);
% 
%         delta = delta / eta;
% 
%     end
% 
%     inner_num_catalyst_screening = inner_num;
%     ritr_Catalyst_Screening = ritr / 2;
% 
%     Err_SVRG(exp_time) = err_nested_SVRG_hist(ceil(3*Itrmax/4),1);
%     Err_SVRG_Screening(exp_time) = err_nested_SVRG_Screening_hist(ceil(3*Itrmax/4),1);
%     Err_SPDHG(exp_time) = err_SPDHG_hist(Itrmax,1);
%     Err_Catalyst(exp_time) = err_nested_Catalyst_hist(ceil(Itrmax^2/ritr_Catalyst),1);
%     Err_Catalyst_Screening(exp_time) = err_nested_Catalyst_Screening_hist(ceil(Itrmax^2/ritr_Catalyst_Screening),1);
% 
% end

% Err = [Err_SVRG;Err_SVRG_Screening;Err_SPDHG;Err_Catalyst;Err_Catalyst_Screening];
% 
% g1 = repmat({'Nested-SVRG'},Exp_Time,1);
% g2 = repmat({'Nested-SVRG-Screening'},Exp_Time,1);
% g3 = repmat({'Nested-SGDM'},Exp_Time,1);
% g4 = repmat({'Nested-Catalyst'},Exp_Time,1);
% g5 = repmat({'Nested-Catalyst-Screening'},Exp_Time,1);
% g = [g1;g2;g3;g4;g5];
% boxplot(Err,g)
% xlabel('Methods')
% ylabel('Relative Error After 1e6 Iterations')

% save('Err_Exp4.mat','Err')

% load('Err_Exp4.mat')
% 
% Exp_Time = 10;
% 
% g1 = repmat({'Nested-SVRG'},Exp_Time,1);
% g2 = repmat({'Nested-SVRG-Screening'},Exp_Time,1);
% g3 = repmat({'Nested-SGDM'},Exp_Time,1);
% g4 = repmat({'Nested-Catalyst'},Exp_Time,1);
% g5 = repmat({'Nested-Catalyst-Screening'},Exp_Time,1);
% g = [g1;g2;g3;g4;g5];
% boxplot(log10(Err),g)
% xlabel('Methods')
% ylabel('lg(Relative Error) After 1e6 Iterations')



% load('Err_Exp4.mat')
% 
% Exp_Time = 10;
% 
% g1 = repmat({'Nested-SVRG'},Exp_Time,1);
% % g2 = repmat({'Nested-SVRG-Screening'},Exp_Time,1);
% g3 = repmat({'SPDHG'},Exp_Time,1);
% g4 = repmat({'Nested-Catalyst'},Exp_Time,1);
% % g5 = repmat({'Nested-Catalyst-Screening'},Exp_Time,1);
% g = [g1;g3;g4];
% I = [1:10,21:30,31:40];
% boxplot(log10(Err(I)),g)
% xlabel('Methods')
% ylabel('lg(Relative Error) After 1e6 Iterations')


load('Err_Exp4.mat')

Exp_Time = 10;

g1 = repmat({'Nested-SVRG'},Exp_Time,1);
% g2 = repmat({'Nested-SVRG-Screening'},Exp_Time,1);
g3 = repmat({'SPDHG'},Exp_Time,1);
g4 = repmat({'Nested-Catalyst'},Exp_Time,1);
% g5 = repmat({'Nested-Catalyst-Screening'},Exp_Time,1);
g = [g1;g3;g4];
I = [1:10,21:30,31:40];
boxplot(log10(Err(I)),g)
xlabel('Methods')
ylabel('lg(Relative Error)')
