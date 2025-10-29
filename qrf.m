    function q = qrf(alpha) % quantile function
    
    mu = 1;
    q = 2.*mu.*sqrt((-log(1-alpha))/pi);

    end 