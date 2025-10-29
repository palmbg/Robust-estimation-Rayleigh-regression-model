
function [coef, mu_hat, resid, pvalues, vcov, stderror, R2, aic, bic] = rrfit_r(x,y,p)

x = x';
y = y';
n_y = size(y);
y1 = ones(n_y(1,1),1);
x1 = [y1 x];
n_x1 = size(x1);
r = n_x1(1,2);

n = n_y(1,1);

ystar = log(y);

reg = regress(ystar,x1);
ini = reg'; % chute inicial

epson = 0.1;


loglik1 = makeloglik(double(x1), double(y));
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'Display', 'none');
[opt21, max] = fminunc(loglik1, double(ini), options);

par_in = opt21;

%while(epson > 0.001)

coef1 = par_in;
eta_hat1 = x1 * coef1';
mu_hat1 = exp(eta_hat1);

F = pr(y, mu_hat1)';

w = zeros(1,n);

for i = 1:n
    
    if(F(1,i) < p)
        
        
        w(1,i) = F(1,i) ./ p;
        
    end
    
    if(p <= F(1,i) && F(1,i) <= 1-p)
        
        
        w(1,i) = 1;
        
    end
    
    if(F(1,i) > 1-p)
        
        
        w(1,i) = (1-F(1,i)) ./ p;
        
    end
    
end

vw = diag(w);

loglik = makeloglikr(double(x1), double(y), double(vw));
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'Display', 'none');
[opt2,fval, ~,~,~,hessian] = fminunc(loglik,double(ini), options);

% par_out = opt2;
% 
% epson1 = sort((abs((par_out - par_in) ./ par_in)));
% 
% epson = epson1(1,length(epson1));
% 
% par_in = par_out;
% 
% end

coef = opt2;
eta_hat = x1 * coef';
mu_hat = exp(eta_hat);
W = diag(((4)/(mu_hat.^2))*(mu_hat.^2));
K = x1' * W * x1;
vcov = inv(chol(K)) * (inv(chol(K)))';

cumulativef = pr(y,mu_hat);

resid = norminv(cumulativef);

stderror = sqrt(diag(vcov));

zstat = abs(coef./stderror');

pvalues = 2*(1 - normcdf(zstat));

aic = 2.*fval+2.*(r);
%
bic = 2.*fval+log(n).*(r);

ini2 = mean(y);

loglik1 = logliknull(double(y));
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', false, 'Display', 'none');
[opt3,fval1] = fminunc(loglik1, double(ini2), options);

mu_hat_n = opt3;

F1 = pr(y, mu_hat_n)';

w1 = zeros(1,n);

for i = 1:n
    
    if(F1(1,i) < p)
        
        
        w1(1,i) = F1(1,i) ./ p;
        
    end
    
    if(p <= F1(1,i) && F1(1,i) <= 1-p)
        
        
        w1(1,i) = 1;
        
    end
    
    if(F1(1,i) > 1-p)
        
        
        w1(1,i) = (1-F1(1,i)) ./ p;
        
    end
    
end

vw1 = diag(w1);


loglik1 = logliknullr(double(y),double(vw1));
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', false, 'Display', 'none');
[~,fval1] = fminunc(loglik1, double(ini2), options);


R2 = 1 - exp(-((2/n)) .* (-fval + fval1));



end


