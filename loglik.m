   function ll = loglik(beta)
   
    eta1 = x1*beta;

    mu = exp(eta1); 
  
    ll = sum(log(pi/2)+log(y)-log(mu^2)-(pi*y^2)/(4*(mu^2)));
    
   end
  