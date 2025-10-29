function ll = makeloglik_star(x1, y)
ll = @loglik;

    function [yy, J] = loglik(beta)
        
        mu = exp(x1*beta');
        
        yy = -sum(sum(log(pi/2)+log(y)-log(mu.^2)-(pi.*y.^2)./(4.*(mu.^2))));
        
        J_star = (x1' * diag(mu) * ((pi.*(y.^2))./(2.*(mu.^3))-(2)./(mu)))';
        
        nm = max(size(beta));
        
        W = diag(((4)/(mu.^2))*(mu.^2));
        K = x1' * W * x1;
        vcov = inv(chol(K)) * (inv(chol(K)))';
        
        
        a = diag ( x1 * vcov * x1');
        
        W1 = diag(-4);
        
        bias_star = eye(nm,nm) * vcov * x1' * W1 * a;
        
        J = -(J_star - (bias_star' * K));
        
    end

end