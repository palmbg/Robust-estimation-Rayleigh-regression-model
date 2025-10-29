function ll = logliknullr(y, vw1)
ll = @loglik;

    function [yy J] = loglik(beta)
        
        mu = beta;
    
        yy =  -sum(vw1*(log(pi/2)+log(y)-log(mu.^2)-(pi.*y.^2)./(4.*(mu.^2))));
        %yy = -(w1*(log(pi/2)+log(y)-log(mu'.^2)-(pi.*y.^2)./(4.*(mu'.^2))));
        %
        J = -sum(sum((vw1 * pi .* y.^2 - 4 .* mu.^2) ./ (2 .* mu.^3)));
        
        
    end

end