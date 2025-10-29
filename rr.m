    function y = rr(mu) % metodo da inversao
        
    n = size(mu);
    u = rand(1,n(1,2));
    y = 2.*mu.*sqrt(-log(1-u)/pi); % metodo da inversao
        
    end 