function yy = Objectivefunction1(x1, y, beta)

        mu = exp(x1*beta');
        
        yy = -sum(sum(log(pi/2)+log(y)-log(mu.^2)-(pi.*y.^2)./(4.*(mu.^2))));
        
  
end


