rng(0)
s = rng;

mnum = 52;

nnum = mnum;

knum = 24;

tnum = 10;

LA = 2;

mu = 1;

H = mu * eye(nnum);

A = zeros(mnum,nnum);

b = zeros(mnum,1);

xst0 = zeros(nnum,1);

yst0 = zeros(mnum,1);

for i = 1 : 2 * knum
    xst0(i) = i;
    yst0(i) = mu / LA * i * ( 4 * knum - i + 1 );
    b(i) = LA / 2;
    A( i, 2 * knum + 1 - i ) = LA / 2;
    if i == 2 * knum
        break
    end
    A( i, 2 * knum - i ) = -LA / 2;
end

for i = 2 * knum + 1 : mnum
    A(i,i) = LA;
end

Atemp = A;

% cond(A) 

xi0 = norm(yst0,"inf");

ratio = [1;2;4;8];

ynrm = zeros(4,1);

for itrnum = 1 : 4

    xi = ratio(itrnum)*xi0;

    epsilon = sqrt(5) * LA * norm(yst0,2) / mu / 16 / ( 2 * tnum + 5 );

    f = @(delta) 2 * mnum * delta * log(xi/mu/delta) - epsilon;

    delta = fsolve(f,1);

    delta = abs(delta);

    delta = delta * 2;

    L = mu + xi * LA ^ 2 / 2 / delta;

    U = eye(nnum);

    V = eye(mnum);

    x0 = [randn(2*knum,1);zeros(nnum-2*knum,1)];

    K = zeros(nnum, 2 * tnum + 2);

    K(:,1) = b;

    A2 = A*A' / norm(A,2) ^ 2;

    for i = 2 : 2 * tnum + 2
        K(:,i) = A2 * K(:,i-1);
        K(:,i) = K(:,i) / norm(K(:,i),2);
    end

    J = K;

    K = A'*K;

    x = x0;
    z = x;

    y0 = zeros(mnum,1);
    y = zeros(mnum,1);

    for i = 1 : mnum
        y0(i) = xi * ( exp((A(i,:)*x-b(i))/delta) - 1 ) /  ( exp((A(i,:)*x-b(i))/delta) + 1 );
    end

    Phi = newtransform(x,K(:,1),K(:,2));
    Psi = newtransform(y0,J(:,1),J(:,2));

    U = U * Phi;
    V = V * Psi;

    Atemp = Psi' * Atemp * Phi;

    Q = L / mu;

    for i = 1 : tnum

        g = newgrad(x,H,xi,delta,Atemp,b);

        ztemp = x - g / L;

        x = (1 + (sqrt(Q)-1)/(sqrt(Q)+1)) * ztemp - (sqrt(Q)-1)/(sqrt(Q)+1) * z;

        z = ztemp;

        Ktemp = U'* K(:, 1 : 2 * i + 2 );

        Jtemp = V'* J(:, 1 : 2 * i + 2 );

        for t = 1 : mnum
            y(t) = xi * ( exp((Atemp(t,:)*x-b(t))/delta) - 1 ) /  ( exp((Atemp(t,:)*x-b(t))/delta) + 1 );
        end

        if isnan(Phi)
        end

        Phi = newtransform(x,Ktemp(:,1: 2*i+1),Ktemp(:,2*i+2));
        Psi = newtransform(y,Jtemp(:,1: 2*i+1),Jtemp(:,2*i+2));

        U = U * Phi;
        V = V * Psi;

        Atemp = Psi' * Atemp * Phi;


    end


    xf = x;
    yf = y;

    x = x0;
    z = x0;
    y = y0;


    for i = 1 : tnum

        g = newgrad(x,H,xi,delta,Atemp,b);

        ztemp = x - g / L;

        x = (1 + (sqrt(Q)-1)/(sqrt(Q)+1)) * ztemp - (sqrt(Q)-1)/(sqrt(Q)+1) * z;

        z = ztemp;

    end

    min_A_val = 2;
    for i = 1 : mnum
        min_A_val = min(min_A_val,norm(Atemp(i,:),2));
    end

    ynew = V' * yst0;
    ynrm(itrnum) = norm(ynew,'inf');

end

plot(1:4,xi0*ratio,1:4,ynrm)
legend('xi','\|y\|_{\infty}')
