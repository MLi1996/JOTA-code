%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 Example of Catalyst/QuickeNing-SVRG                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Part 1:  load and prepare data                                        %
%   Part 2:  apply SVRG and Catalyst-SVRG                                 %
%   Part 3:  plot comparison                                              %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%   Loading training and testing data    %%%%%%%%%%%%%%%%%%%
%                                                                         %
% DETAILS SEE load_data.m                                                 %
%                                                                         %
% TO USE YOUR OWN DATASET:                                                %
% a) Place the training data in: v1/data/NAME_train.mat                   %
% b) Set the variables names as                                           %
%       Xtrain : p*n matrix (p : dimension of feature, n: size of dataset)%
%       Ytrain : n*1 matrix                                               %
% c) (OPTIONAL) Place the testing data in: v1/data/NAME_test.mat          %
% d) Set the variables names as                                           %
%       Xtest : p*ntest matrix (ntest: size of testing data)              %
%       Ytest : ntest*1 matrix                                            %
%                                                                         %
% IMPORTANT MESSAGE:                                                      %
%       1) X, Y must be arrays of doubles                                 %
%       2) NORMALIZE YOUR DATASET                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dataset = 'covtype';     % NAME of dataset
% [Xtrain,Ytrain,Xtest,Ytest] = load_data(dataset);


load("mushrooms_matrix.mat");
load("mushrooms_vector.mat");

m = 8124;
n = 111;

A_org = full(instance_matrix);

A_org = normalize(A_org);
A_org = [A_org(:,1:77),A_org(:,79:112)];
A = A_org;
Xtrain0 = A';
Ytrain0 = label_vector;

for i = 1 : m
    if Ytrain0(i) == -1
        Xtrain0(:,i) = - Xtrain0(:,i);
        Ytrain0(i) = 1;
    end
end

Xtest = [];
Ytest = [];


% % Trivial Test
% Xtrain = [1,-1];
% Ytrain = [1;-1];
% 
% Xtest = [];
% Ytest = [];

%%%%%%%%%%%%%%%%%%%%   Specify the LOSS FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Available loss functions:                                               %
% a) 'logi' : logistic loss + l2 regularization                           %
%       l2 parameter: mu/n                                                %
% b) 'elasticnet': square loss + l2 regularization + l1 regularization    %
%       l2 parameter: mu/n                                                %
%       l1 parameter: lambda/n                                            %
% c) 'lasso': square loss + l1 regularization                             %
%       l1 parameter: lambda/n                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 'logi';     % Models: 'logi', 'elasticnet' or 'lasso'


xi = 0.5;
delta0 = 0.02;


b = -ones(m,1);

for i = 1 : m
    A(i,:) = - A(i,:) * label_vector(i);
end

cvx_begin
    variable x(n,1);
    dual variable lambda;
    minimize 1 / 2 * (x'* x);
    subject to
        lambda: A * x <= b;
cvx_end

x_cvx = x;

delta = delta0;

x0 = zeros(n,1);

x = x0;

num = 0;

nummax = 1250;

err_list = [];

nb = 100;

while num <= nummax

    Xtrain = Xtrain0;
    Ytrain = Ytrain0 / delta;
    mu = 1 / xi / delta;
    default_catalyst = 1;   % Use default parameter setting
    algo = 'svrg';          % Algo in the inner loop of Catalyst
    % For the current version, only svrg is available

    [param] = param_catalyst(delta,Xtrain, model, mu, lambda, algo, default_catalyst);

    nb_it = nb;            % Nb of epochs

    [w_catalyst,it_catalyst,train_loss_list_catalyst,dualgaplist_catalyst, nnzlist_catalyst, err_temp]...
        = catalyst_withouttest(x,x_cvx,Xtrain,Ytrain,param, nb_it);
    
    x = w_catalyst;
    
    err_list = [err_list,err_temp(2:length(err_temp))];
    

    num = num + nb;

    delta = delta / 2;
    nb = nb * 2;

end




semilogy(m*(1:length(err_list)), err_list);
xlabel('Number of Gradient Evaluations');
ylabel('Relative Error');


save("Catalyst_Err_List.mat","err_list");

















% fun = @(w) -1 / (1+exp(w-1)) - 1 / (1+exp(w+1)) + w ; 
% w_mine = fsolve(fun,0);



% val1 = norm(w_catalyst-x_cvx,2) / norm(x_cvx);
% w = w_catalyst;
% grad_res = w;
% 
% for i = 1 : m
%     grad_res = grad_res + xi * exp((A(i,:)*w-b(i))/delta) / (1+exp((A(i,:)*w-b(i))/delta)) * (A(i,:))'; 
% end
% 

% v1 = label_vector / delta ;
% v1 = v1';
% 
% A1 = Xtrain0 * diag(label_vector) /delta;
% cvx_solver sedumi
% cvx_expert true
% cvx_begin
%     variables x(n)
%     minimize( 1/2 *x'*x+ xi * delta*sum(log_sum_exp([zeros(1,m); v1-x'*A1])))
% cvx_end
% 
% x_pen = x;
% val2 = norm(w_catalyst-x_pen,2) / norm(x_pen);



