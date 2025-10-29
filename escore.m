   function Ubeta = escore(beta)
  
    eta1 = x1*beta;

    mu = exp(eta1);
    
    T = diag(mu);
    
    v = (pi.*(y.*y))/(2.*(mu.*mu.*mu))-(2)/(mu);
    
    Ubeta  = t(x1);
   
   end