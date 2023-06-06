function [xk,it,train_loss_list,dualgaplist,nnzlist,err_list] = catalyst_withouttest(x0,xst,Xtrain,Ytrain,param, nb_it)

%   This function run catalyst to optimize the given objective function
%   Outputs:
%       xk              : the last iterate
%       it              : list of nb of epochs passed 
%       train_loss_list : list of values of f(xk) 
%       dualgaplist     : list of Fenchel duality gap at each xk
%       train_acc_list  : list of training accuracy at each xk
%       test_acc_list   : list of testing accuracy at each xk
%       test_loss_list  : list of loss functions on test set at each xk 
%       nnzlist         : list of number of non zero components of xk

ntrain=size(Xtrain,2);
p=size(Xtrain,1);


%%%% Set parameters for Test evaluation
param_test = param;
param_test.mu = 0;
param_test.lambda = 0;

%%%% Initialization
xk = x0;
train_loss= compute_loss(xk,Ytrain,Xtrain,param);
train_loss_list = [train_loss];

err_list = [norm(xk-xst,2)/norm(xst,2)];

it = [0];

%%%% Evaluate Train and Test accuracy
Ytrain_pred = linear_prediction(Xtrain,xk);
train_acc = sum(Ytrain_pred == Ytrain)/ntrain;



%%%% Evaluate Duality gap
dualgap = compute_dualgap(xk,Ytrain,Xtrain,param); % Fenchel duality gap
dualgap_init = min(dualgap,train_loss);
fprintf('Inital duality gap :%g \n',dualgap_init);

dualgaplist = [dualgap];
nnzlist = [nnz(xk)];
total_it =0;
count_it =0;

q = param.mu/(param.mu+param.kappa);
rho = 0.9*sqrt(q);
tau = 1- rho;
gamma_power = 0.1;
yk = xk;
yk_past = yk;

if param.mu == 0
   alpha = 1;  
else
   alpha = sqrt(q);
end

while total_it < nb_it
    count_it = count_it+1;
    if strcmp(param.stop_criterion,'absolute')
        if param.mu > 0
            param.epsilon = max(2*dualgap_init*tau^count_it/9,10^(-15));
        else
            param.epsilon = max(2*dualgap_init/(9*(count_it+1)^(4+gamma_power)),10^(-15));
        end
    elseif strcmp(param.stop_criterion,'relative')
        if param.mu > 0 
            param.delta = max( sqrt(q)/(2-sqrt(q)),10^(-10));
        else 
            param.delta = max(1/((count_it+1)^2),10^(-10)); 
        end
    end
    
    if strcmp(param.algo,'svrg')
        [g, F,xk_new, nb_grad] = approxgrad(yk,Ytrain,Xtrain,param,xk,yk_past);
    end
    
    total_it = total_it + nb_grad; % Total nb of epochs passed
    alpha_new = (q-alpha^2 + sqrt((q-alpha^2)^2+4*alpha^2))/2;
    beta = (alpha*(1-alpha))/(alpha^2+alpha_new);
    yk_past = yk;
    yk = xk_new + beta*(xk_new - xk);
    %%%% Renew variables
    xk = xk_new;
    alpha = alpha_new;


    %%%% Save values
    train_loss = compute_loss(xk,Ytrain,Xtrain,param);
    if param.mu > 0 || param.lambda >0
        dualgap = compute_dualgap(xk,Ytrain,Xtrain,param);
        fprintf('Iter %d, loss: %g, dualgap: %g \n',total_it, train_loss,dualgap);
        dualgaplist= [dualgaplist, dualgap];
    else
        fprintf('Iter %d, Train loss: %g \n',total_it, train_loss);
    end
    
    %%%% Train accuracy 
    Ytrain_pred = linear_prediction(Xtrain,xk);
    train_acc = sum(Ytrain_pred == Ytrain)/ntrain;
    train_loss_list = [train_loss_list, train_loss]; 
    err_list = [err_list,norm(xk-xst,2)/norm(xst,2)];
    it = [it, total_it];
    nnzlist = [nnzlist, nnz(xk)];

    


end


end