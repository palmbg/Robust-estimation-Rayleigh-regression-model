n = 150; %tamanho amostral
R = 5000;
%
% %parametros
beta0 = 0.5;
beta1 = 0.15;

%
p = 0.001;
par = [beta0,beta1];
%
    % %covariaveis
    x1 = rand(1,n);
    %
    eta = beta0+beta1*x1; % preditor linear
    mu = exp(eta); % media
    %
    
    y = rr(mu); % variavel dependente
    x = [x1];
    
    y1 = y;
    
    y1(1,1:r) = y1(1,1:r) + 10;
    %
    
    rrfit_r(x,y1,p)
