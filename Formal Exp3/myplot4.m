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