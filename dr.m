   function d = dr(x,mu) % density function
  
   %mu = 1;
   d = pi.*x ./(2.*mu.^2).*exp(-(pi.*(x.^2)) ./(4.*mu.^2));
  
   end 