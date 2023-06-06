rng(1)

l = 100;
m = 100;
n = 100;

T = 1;

Mul = [0.36; 0.6; 1];
Eta = [2; 4; 8];
Alpha = [0.5; 0.8; 0.9];

err_SASC = zeros(T,length(Mul),length(Eta));
err_nested_SGD = zeros(T,length(Mul),length(Eta));
err_nested_SGDM = zeros(T,length(Mul),length(Eta),length(Alpha));

for t = 1 : T

    A = randn(m,n);
    Phi = randn(l,n);
    y = randn(l,1);
    b = abs(randn(m,1));

    w = 0.1;

    xi = 1;

    x0 = zeros(n,1);

    for i = 1 : m
        A(i,:) = A(i,:) / norm(A(i,:));
    end

    delta0 = 0.05;

    cvx_begin
        variable x(n,1);
        dual variable lambda;
        minimize (1/2/l * sum_square(Phi*x-y) + w / 2 * x'* x);
        subject to
            lambda: A * x <= b;
    cvx_end

    x_cvx = x;
    lambda_cvx = lambda;

    Itrmax = 1e7 + 1;

    sample = unidrnd(m,Itrmax,1);

    L = w + norm(Phi,2) ^ 2 / l;

    mu = w;

    % SASC-SGD
    for mulnum = 1 : length(Mul)
        for etanum = 1 : length(Eta)
            x = x0;
            mul = Mul(mulnum);
            eta = Eta(etanum);
            itr = 0;

            inner_num = 0;

            s = 1 / 8 / L;

            while inner_num >= 0

                inner_num = inner_num + 1;

                beta = 4 * s;

                M = ceil(mul * eta / mu / s );

                for i = 1 : M
                    itr = itr + 1;
                    j = sample(itr);

                    g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + (A(j,:))'*max(A(j,:)*x-b(j),0) / beta;

                    x = x - s * g;

                    if itr >= Itrmax
                        break
                    end

                end

                if itr >= Itrmax
                    break
                end

                s = s / eta;

            end
            err_SASC(mulnum,etanum) = norm(x-x_cvx,2)/norm(x_cvx,2);
        end
    end


    % Nested-SGD

    for mulnum = 1 : length(Mul)
        for etanum = 1 : length(Eta)
            x = x0;
            mul = Mul(mulnum);
            eta = Eta(etanum);
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

                    g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

                    x = x - s * g;

                    if itr >= Itrmax
                        break
                    end

                end

                if itr >= Itrmax
                    break
                end

                delta = delta / eta;

            end
            err_nested_SGD(mulnum,etanum) = norm(x-x_cvx,2)/norm(x_cvx,2);
        end
    end

    

    % Nested-SGDM
    for mulnum = 1 : length(Mul)
        for etanum = 1 : length(Eta)
            for alphanum = 1:length(Alpha)
                x = x0;
                mul = Mul(mulnum);
                eta = Eta(etanum);
                alpha = Alpha(alphanum);
                v = zeros(n,1);
                itr = 0;
                inner_num = 0;

                delta = delta0;

                while inner_num >= 0

                    inner_num = inner_num + 1;

                    s = 1 / ( L + mu + m * xi / 4 / delta );

                    M = ceil( mul * ( 2 * log(2*eta-1) + log( L / mu + m * xi / 4 / mu / delta ) ) * sqrt( m * ( L / mu + m * xi / 4 / mu / delta ) ) );

                    v = zeros(n,1);

                    for i = 1 : M
                        itr = itr + 1;
                        j = sample(itr);

                        g = ( Phi(j,:) )' * ( Phi(j,:) * x - y(j) ) + w * x + m * xi * ( A(j,:) )' * expd((A(j,:)*x-b(j))/delta);

                        v = alpha * v - s * g;

                        x = x + v;

                        if itr >= Itrmax
                            break
                        end

                    end

                    if itr >= Itrmax
                        break
                    end

                    delta = delta / eta;
                    
                end
                err_nested_SGDM(mulnum,etanum,alphanum) = norm(x-x_cvx,2)/norm(x_cvx,2);
            end
        end
    end

    


end
       