% rng(5)
% 
% num_C = 4;
% 
% n = 500;
% l = 100;
% m = 100;
% 
% 
% W = randn(num_C,n);
% 
% for i = 1:num_C
%     W(i,:) = W(i,:) / norm(W(i,:));
% end
% 
% X10 = randn(n,l);
% Z1 = zeros(1,l);
% Y1 = zeros(num_C,l);
% 
% for i = 1:l
%     X10(:,i) = X10(:,i) / norm(X10(:,i));
%     [~, idxtemp] = max(W*X10(:,i));
%     Z1(i) = idxtemp;
%     Y1(idxtemp,i) = 1;
% end
% 
% X2 = randn(n,m);
% Z2 = zeros(1,m);
% Y2 = zeros(num_C,m);
% 
% mytest = zeros(1,m);
% 
% for i = 1:m
%     X2(:,i) = X2(:,i) / norm(X2(:,i));
%     temp = W*X2(:,i);
%     temp = sort(temp);
%     mytest(i) = temp(4)-temp(3);
%     [~, idxtemp] = max(W*X2(:,i));
%     Z2(i) = idxtemp;
%     Y2(idxtemp,i) = 1;
% end
% 
% X11 = randn(n,l);
% 
% 
% for i = 1:l
%     X11(:,i) = X11(:,i) / norm(X11(:,i));
% end
% 
% X1 = X10 + 0.1 * X11;
% 
% 
% mu = 0.05;
% 
% 
% A2 = zeros( (num_C-1) * m , num_C * n );
% 
% for i = 1 : m
%     index = ( num_C - 1 ) * ( i - 1 );
%     pos = Z2(i);
%     pos_non = setdiff(1:num_C,pos); 
%     for j = 1 : num_C - 1
%         A2(index + j, ( (pos-1) * n + 1 ) : (pos * n) ) = (X2(:,i))';
%         posn = pos_non(j);
%         A2(index + j, ( (posn-1) * n + 1 ) : (posn * n) ) = -(X2(:,i))';
%     end
% end
% 
% 
% cvx_begin
%     variable w(n,num_C);
%     dual variable lambda;
%     minimize (-(vec(w*Y1))'*(vec(X1)) + sum(log_sum_exp((w')*X1)) + mu / 2 * sum_square(vec(w)));
%     subject to
%         lambda: A2 * vec(w) >= 1;
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
% mul_SVRG = 1;
% eta_SVRG = 2;
% 
% mul_SVRG_screening = 1;
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
% 
% S = kron((eye(num_C)-ones(num_C,num_C)/num_C), X1*X1') / 2;
% 
% L = max(eig(S)) + mu;
% 
% mu = mu;
% 
% xi = 1;
% delta0 = 1;
% 
% 
% w0 = zeros(n,num_C);
% Itrmax1 = Itrmax;
% Itrmax2 = Itrmax;
%        
% sample2 = unidrnd(m,Itrmax2,1);
% 
% SVRG_ratio = 3;
% Catalyst_Budget = 1;
% 
% 
% 
% % Nested-Catalyst_Screening
% 
% err_nested_Catalyst_hist = zeros(Itrmax2,1);
% 
% w = w0;
% mul = 1.5;
% eta = 2;
% N = SVRG_ratio;
% 
% T_B = Catalyst_Budget;
% 
% ritr = 0;
% 
% mnewtemp = (num_C-1) * m; 
% 
% Itemp = 1:mnewtemp;
% 
% itr = 0;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     wy = w;
% 
%     wold = w;
% 
%     inner_num = inner_num + 1;
% 
%     Ltemp = L + mnewtemp / 3 * xi / 4 / delta;
% 
%     kappa = ( L - mu ) / m ;
% 
%     q = mu / ( mu + kappa );
% 
%     alpha = sqrt(q);
% 
%     beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );
% 
%     s = 2 / ( Ltemp + mu + kappa );
% 
%     M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + mnewtemp / 3 * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + (mnewtemp / 3) ^2 * xi / 4 / mu / delta  ) );
% 
%     M = ceil( M / ( N + 1 ) / mnewtemp / T_B) * ( N + 1 ) * mnewtemp * T_B;
% 
% 
%     sample_temp = unidrnd(mnewtemp,3*M,1);
% 
%     i = 0;
% 
%     obj_i = N * l;
% 
%     cons_i = ceil( N * mnewtemp / 3 );
%     cons_j = 0;
%     
%     t_b = 0;
% 
%     while i < 6 * M % 3 = num_C-1
%         if obj_i >= N * l
% 
%             y = w;
%             v1 = zeros(n,num_C);
% 
%             for k = 1 : m
%                 g = zeros(n,num_C);
% 
%                 g(:,Z1(k)) = g(:,Z1(k)) - X1(:,k);
% 
%                 myexp1 = exp(w'*X1(:,k));
% 
%                 for kk = 1 : num_C
%                     g(:,kk) = g(:,kk) + myexp1(kk) / sum(myexp1) * X1(:,k);
%                 end
%                 v1 = v1 + g;
%             end
% 
%             ritr = ritr + 3 * m;
%             i = i + 3 * m;
%             obj_i = 0;
% 
%         end
% 
%         if cons_i >= N * mnewtemp/3
% 
%             z = w;
%             v2 = zeros(n,num_C);
% 
%             for ktemp = 1 : mnewtemp
% 
%                 g = zeros(n,num_C);
%                 index_temp = Itemp(ktemp);
%                 
%                 k = ceil(index_temp/3);
% 
%                 kk = index_temp - ( k - 1 ) * 3;
% 
%                 pos = Z2(k);
% 
%                 pos_non = setdiff(1:num_C,pos);
% 
%                 val = expd( ((X2(:,k))' * ( w(:,pos_non(kk)) - w(:,pos) )+1) / delta ) * xi;
%                 g(:,pos_non(kk)) = g(:,pos_non(kk)) + val * X2(:,k);
%                 g(:,pos) = g(:,pos) - val * X2(:,k);
% 
%                 v2 = v2 + g;
% 
%             end
% 
%             ritr = ritr + mnewtemp;
%             i = i + mnewtemp;
%             cons_i = 0;
% 
%         end
% 
%         
%         itr = itr + 1;
%         i = i + 6;
%         ritr = ritr + 6;
%         obj_i = obj_i + 1;
%         cons_i = cons_i + 1;
%         cons_j = cons_j + 1;
%         J_temp = sample_temp((3*(cons_j-1)+1):(3*cons_j));
%         j = sample2(itr);
%         
%         gx_obj = zeros(n,num_C);
%         gy_obj = zeros(n,num_C);
%         
%         gx_obj(:,Z1(j)) = gx_obj(:,Z1(j)) - l * X1(:,j);
% 
%         myexp1 = exp(w'*X1(:,j));
% 
%         for k = 1 : num_C
%             gx_obj(:,k) = gx_obj(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j); 
%         end
%         
%         gy_obj(:,Z1(j)) = gy_obj(:,Z1(j)) - l * X1(:,j);
% 
%         myexp1 = exp(y'*X1(:,j));
% 
%         for k = 1 : num_C
%             gy_obj(:,k) = gy_obj(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j); 
%         end
% 
%         gdiff_obj = gx_obj - gy_obj;
% 
%         gx_cons = zeros(n,num_C);
% 
%         gz_cons = zeros(n,num_C);
% 
%         for numtemp = 1 : 3
%             index_temp = Itemp(J_temp(numtemp));
% 
%             k = ceil(index_temp/3);
% 
%             kk = index_temp - ( k - 1 ) * 3;
% 
%             pos = Z2(k);
% 
%             pos_non = setdiff(1:num_C,pos);
% 
%             val = expd( ((X2(:,k))' * ( w(:,pos_non(kk)) - w(:,pos) )+1) / delta ) * xi * mnewtemp / 3 ;
%             gx_cons(:,pos_non(kk)) = gx_cons(:,pos_non(kk)) + val * X2(:,k);
%             gx_cons(:,pos) = gx_cons(:,pos) - val * X2(:,k);
% 
%             val = expd( ((X2(:,k))' * ( z(:,pos_non(kk)) - z(:,pos) )+1) / delta ) * xi * mnewtemp / 3 ;
%             gz_cons(:,pos_non(kk)) = gz_cons(:,pos_non(kk)) + val * X2(:,k);
%             gz_cons(:,pos) = gz_cons(:,pos) - val * X2(:,k);
%         end
% 
%         gdiff_cons = gx_cons - gz_cons;
% 
%         w = ( w - s * gdiff_obj - s * v1 - s * gdiff_cons - s * v2 ) / ( 1 + s * mu + s * kappa ) + wy * s * kappa / ( 1 + s * mu + s * kappa );
% 
%         err_nested_Catalyst_hist(itr) = norm(w-w_cvx,'fro') / norm(w_cvx,'fro');
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%         if cons_i >= N * mnewtemp/3
%             t_b = t_b + 1;
%         end
%         
%         if t_b == T_B
%             wy = w + beta * (w - wold);
%             wold = w;
%             t_b = 0;
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
% %     Itemp_new = [];
% % 
% %     for i = 1 : mnewtemp
% %         ii = Itemp(i);
% %         k = ceil(ii/3);
% % 
% %         kk = ii - ( k - 1 ) * 3;
% % 
% %         pos = Z2(k);
% % 
% %         pos_non = setdiff(1:num_C,pos);
% % 
% %         val = (X2(:,k))' * ( w(:,pos_non(kk)) - w(:,pos) )+1;
% % 
% %         if val >= -2 * sqrt(mnewtemp) * delta
% %             Itemp_new = [Itemp_new,ii];
% %         end
% %     end
% % 
% %     Itemp = Itemp_new;
% %     mnewtemp = length(Itemp);
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_Catalyst = inner_num;
% ritr_num_Catalyst = ritr/6;
% 
% 
% 
% 
% 
% 
% 
% % Nested-Catalyst_Screening
% 
% err_nested_Catalyst_Screening_hist = zeros(Itrmax2,1);
% 
% w = w0;
% mul = 1.5;
% eta = 2;
% N = SVRG_ratio;
% 
% T_B = Catalyst_Budget;
% 
% ritr = 0;
% 
% mnewtemp = (num_C-1) * m; 
% 
% Itemp = 1:mnewtemp;
% 
% itr = 0;
% 
% inner_num = 0;
% 
% delta = delta0;
% 
% while inner_num >= 0
% 
%     wy = w;
% 
%     wold = w;
% 
%     inner_num = inner_num + 1;
% 
%     Ltemp = L + mnewtemp / 3 * xi / 4 / delta;
% 
%     kappa = ( L - mu ) / m ;
% 
%     q = mu / ( mu + kappa );
% 
%     alpha = sqrt(q);
% 
%     beta = alpha * ( 1 - alpha ) / ( alpha^2 + alpha );
% 
%     s = 2 / ( Ltemp + mu + kappa );
% 
%     M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + mnewtemp / 3 * xi / 4 / mu / delta ) ) * sqrt(  l * L / mu + (mnewtemp / 3) ^2 * xi / 4 / mu / delta  ) );
% 
%     M = ceil( M / ( N + 1 ) / mnewtemp / T_B) * ( N + 1 ) * mnewtemp * T_B;
% 
% 
%     sample_temp = unidrnd(mnewtemp,3*M,1);
% 
%     i = 0;
% 
%     obj_i = N * l;
% 
%     cons_i = ceil( N * mnewtemp / 3 );
%     cons_j = 0;
%     
%     t_b = 0;
% 
%     while i < 6 * M % 3 = num_C-1
%         if obj_i >= N * l
% 
%             y = w;
%             v1 = zeros(n,num_C);
% 
%             for k = 1 : m
%                 g = zeros(n,num_C);
% 
%                 g(:,Z1(k)) = g(:,Z1(k)) - X1(:,k);
% 
%                 myexp1 = exp(w'*X1(:,k));
% 
%                 for kk = 1 : num_C
%                     g(:,kk) = g(:,kk) + myexp1(kk) / sum(myexp1) * X1(:,k);
%                 end
%                 v1 = v1 + g;
%             end
% 
%             ritr = ritr + 3 * m;
%             i = i + 3 * m;
%             obj_i = 0;
% 
%         end
% 
%         if cons_i >= N * mnewtemp/3
% 
%             z = w;
%             v2 = zeros(n,num_C);
% 
%             for ktemp = 1 : mnewtemp
% 
%                 g = zeros(n,num_C);
%                 index_temp = Itemp(ktemp);
%                 
%                 k = ceil(index_temp/3);
% 
%                 kk = index_temp - ( k - 1 ) * 3;
% 
%                 pos = Z2(k);
% 
%                 pos_non = setdiff(1:num_C,pos);
% 
%                 val = expd( ((X2(:,k))' * ( w(:,pos_non(kk)) - w(:,pos) )+1) / delta ) * xi;
%                 g(:,pos_non(kk)) = g(:,pos_non(kk)) + val * X2(:,k);
%                 g(:,pos) = g(:,pos) - val * X2(:,k);
% 
%                 v2 = v2 + g;
% 
%             end
% 
%             ritr = ritr + mnewtemp;
%             i = i + mnewtemp;
%             cons_i = 0;
% 
%         end
% 
%         
%         itr = itr + 1;
%         i = i + 6;
%         ritr = ritr + 6;
%         obj_i = obj_i + 1;
%         cons_i = cons_i + 1;
%         cons_j = cons_j + 1;
%         J_temp = sample_temp((3*(cons_j-1)+1):(3*cons_j));
%         j = sample2(itr);
%         
%         gx_obj = zeros(n,num_C);
%         gy_obj = zeros(n,num_C);
%         
%         gx_obj(:,Z1(j)) = gx_obj(:,Z1(j)) - l * X1(:,j);
% 
%         myexp1 = exp(w'*X1(:,j));
% 
%         for k = 1 : num_C
%             gx_obj(:,k) = gx_obj(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j); 
%         end
%         
%         gy_obj(:,Z1(j)) = gy_obj(:,Z1(j)) - l * X1(:,j);
% 
%         myexp1 = exp(y'*X1(:,j));
% 
%         for k = 1 : num_C
%             gy_obj(:,k) = gy_obj(:,k) + l * myexp1(k) / sum(myexp1) * X1(:,j); 
%         end
% 
%         gdiff_obj = gx_obj - gy_obj;
% 
%         gx_cons = zeros(n,num_C);
% 
%         gz_cons = zeros(n,num_C);
% 
%         for numtemp = 1 : 3
%             index_temp = Itemp(J_temp(numtemp));
% 
%             k = ceil(index_temp/3);
% 
%             kk = index_temp - ( k - 1 ) * 3;
% 
%             pos = Z2(k);
% 
%             pos_non = setdiff(1:num_C,pos);
% 
%             val = expd( ((X2(:,k))' * ( w(:,pos_non(kk)) - w(:,pos) )+1) / delta ) * xi * mnewtemp / 3 ;
%             gx_cons(:,pos_non(kk)) = gx_cons(:,pos_non(kk)) + val * X2(:,k);
%             gx_cons(:,pos) = gx_cons(:,pos) - val * X2(:,k);
% 
%             val = expd( ((X2(:,k))' * ( z(:,pos_non(kk)) - z(:,pos) )+1) / delta ) * xi * mnewtemp / 3 ;
%             gz_cons(:,pos_non(kk)) = gz_cons(:,pos_non(kk)) + val * X2(:,k);
%             gz_cons(:,pos) = gz_cons(:,pos) - val * X2(:,k);
%         end
% 
%         gdiff_cons = gx_cons - gz_cons;
% 
%         w = ( w - s * gdiff_obj - s * v1 - s * gdiff_cons - s * v2 ) / ( 1 + s * mu + s * kappa ) + wy * s * kappa / ( 1 + s * mu + s * kappa );
% 
%         err_nested_Catalyst_Screening_hist(itr) = norm(w-w_cvx,'fro') / norm(w_cvx,'fro');
% 
%         if itr >= Itrmax2
%             break
%         end
% 
%         if cons_i >= N * mnewtemp/3
%             t_b = t_b + 1;
%         end
%         
%         if t_b == T_B
%             wy = w + beta * (w - wold);
%             wold = w;
%             t_b = 0;
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
%     for i = 1 : mnewtemp
%         ii = Itemp(i);
%         k = ceil(ii/3);
% 
%         kk = ii - ( k - 1 ) * 3;
% 
%         pos = Z2(k);
% 
%         pos_non = setdiff(1:num_C,pos);
% 
%         val = (X2(:,k))' * ( w(:,pos_non(kk)) - w(:,pos) )+1;
% 
%         if val >= -2 * sqrt(mnewtemp) * delta
%             Itemp_new = [Itemp_new,ii];
%         end
%     end
% 
%     Itemp = Itemp_new;
%     mnewtemp = length(Itemp);
% 
%     delta = delta / eta;
% 
% end
% 
% inner_num_Catalyst_Screening = inner_num;
% ritr_num_Catalyst_Screening = ritr/6;
% 
% load('Exp3_SVRG_Screening.mat')
% semilogy(1.333333*(1:100:Itrmax2),err_nested_SVRG_Screening_hist(1:100:Itrmax2),...
%    ritr_num_Catalyst / Itrmax2 * (1:100:Itrmax2), err_nested_Catalyst_hist(1:100:Itrmax2),...
%    ritr_num_Catalyst_Screening / Itrmax2 * (1:100:Itrmax2), err_nested_Catalyst_Screening_hist(1:100:Itrmax2));
% legend('Nested-SVRG-Screening','Nested-Catalyst','Nested-Catalyst-Screening');
% xlabel('Number of Stochastic Evaluations');
% ylabel('Relative Error');

log_Err = zeros(3*m,1);

log_Err_Itemp = zeros(mnewtemp,1);

i = 0;
for k = 1 : m
    pos = Z2(k);

    pos_non = setdiff(1:num_C,pos);
    
    for kk = 1 : num_C-1    
        i = i + 1;
        Err_temp = (X2(:,k))' * ( -w_cvx(:,pos_non(kk)) + w_cvx(:,pos) )-1;
        if Err_temp <= 1e-10
            log_Err(i) = -10;
        else
            log_Err(i) = log10(Err_temp);
        end
    end


end

for i = 1 : mnewtemp
    ii = Itemp(i);
    k = ceil(ii/3);

    kk = ii - ( k - 1 ) * 3;

    pos = Z2(k);

    pos_non = setdiff(1:num_C,pos);

    Err_temp = (X2(:,k))' * ( -w_cvx(:,pos_non(kk)) + w_cvx(:,pos) )-1;
    if Err_temp <= 1e-10
        log_Err_Itemp(i) = -10;
    else
        log_Err_Itemp(i) = log10(Err_temp);
    end
end


edges = -10: 1: 2;
h1 = histogram(log_Err,edges);
hold on
h2 = histogram(log_Err_Itemp,edges);
legend('Original','Remaining')
xlabel('lg(Slackness)');
ylabel('Number of Constraints')
