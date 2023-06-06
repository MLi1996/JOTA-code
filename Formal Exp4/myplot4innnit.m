log_Err = zeros(m,1);

% log_Err_Itemp = zeros(mtemp,1);

for k = 1 : m

    Err_temp = (A2(k,:)*w_cvx-1);
    if Err_temp <= 1e-10
        log_Err(k) = -10;
    else
        log_Err(k) = log10(Err_temp);
    end

end

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


edges = -10: 1: 2;
h1 = histogram(log_Err,edges);
xlabel('lg(Slackness)');
ylabel('Number of Constraints')
