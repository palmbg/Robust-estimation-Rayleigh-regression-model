function ll = gradfunction(x1, y)
ll = @escore;

    function [J] = loglik(beta)
        
        mu = exp(x1*beta');
        
        J = (x1' * diag(mu) * ((pi.*(y.^2))./(2.*(mu.^3))-(2)./(mu)))';
    end

 
end