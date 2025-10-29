function ll = Objectivefunction(x1, y,beta)
ll = @loglik;

    function [yy] = loglik(beta)
        
        mu = exp(x1*beta');
        
        yy = -sum(log(pi/2)+log(y)-log(exp((x1*beta')).^2)-(pi.*y.^2) ./ (4.*(exp((x1*beta')).^2)));
        
    end

 
end