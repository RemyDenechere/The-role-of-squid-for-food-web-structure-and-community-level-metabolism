function [R_star , r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f, param, i, j)
% -------------------------------------------------------------------------
% give an approximation of the population growth rate of a population give
% based on the growth coefficient A. 
% -------------------------------------------------------------------------
% This function has been modified from: 
% Title : Deriving population scaling rules from individual-level metabolism and life history traits
% Auteurs : Rémy Denéchère, P Daniël van Denderen, Ken H Andersen
% Date de publication : 2022/4/1
% Revue : The American Naturalist
% Volume : 199
% -------------------------------------------------------------------------
n = 3/4; % Metabolic exponent for fish, the Metabolic exponent for squid is
         % latter set at 2/3 as the data suggest! 

if i == 2 % case 1: Teleost
    %    value         Lower limit    Upper limit 
    A0 = [exp(1.29153), 2.2481,       6.8016];
    b = 0; % slope of the regression 
    ratio = exp(7.299)*Winf.^(0.9085); % M/M_0 
    minWinf = 1; maxWinf = 684000; % g minimum observed teleost
    % conf90_inf = [exp(1.149), 0.02215]; conf90_sup = [exp(1.434), 0.06448]; 
    % 90% confidence is given by the lowest/highest value  A_0 and b (in this
    % order). 
    epsR =  0.03; 
    
% case 2: Elastobranch
elseif i == 3
    A0 = [exp(1.70586), 3.6118, 17.0316]; % former values:  exp(0.5157), exp(2.896)] ; 
    b = [0, -0.09984, 0.1479]; % slope of A relationship
    ratio = 362.5060; % former value : 6.1119*10^3;
    minWinf = 100; maxWinf = 4.8799*10^(5); % g
    % conf90_inf = [exp(0.5157), -0.09984]; conf90_sup = [exp(2.896), 0.1479]; 
    epsR = param.epsR_elasm;    
  
% case 3: Bivalves: 
elseif i == 4
    A0 = [exp(-0.01941), 0.4253, 4.4035]; %  2.2734,11.3766 ]; % former values :exp(-0.05598), exp(0.01716) ;
    b = [0.2149, 0.1888, 0.241] ;
    ratio = exp(14.41)*Winf.^(1.07); % former value: Winf/(5.3*10^(-7));
    minWinf = 6.339*10^(-4); maxWinf = 598.2662; %g 
    % conf90_inf = [exp(-0.05598), 0.1888]; conf90_sup = [exp(0.01716), 0.241]; 
    epsR =  0.03; 
    
% case 4: Copepods active
elseif i == 5  
    A0 = [exp(-2.043)*10^(3/2)/0.6,  exp(1.781), exp(2.064)];
    b = [0]; 
    ratio = 176.6654;  % 439.9322;
    epsR =  0.23; % this is from camila's paper. 
    minWinf = 1.40e-05; maxWinf = 0.0135;
   
% case 5 : Copepods passive
elseif i == 6
    A0 = [exp(0.3652),  exp(0.1141), exp(0.6163)]; % ken's method: 30.3762  125.9025
    b = [0]; 
    ratio = 176.6654;  % 439.9322; % I assume the same ratio for active and 
                                   % passive copepods, i.e., they have the same size 
    epsR =  0.23; 
    minWinf = 3.220e-07; maxWinf = 9.5922e-06;
    
    
% case 6 : unknown organism
elseif i == 1
    A0 = [exp(1.70586) exp(1.70586) exp(1.70586)]; % Based on bivalves
    b = [0.2149]; % Based on bivalves
    ratio = 362.5; % Based on Elasmobranchs  
    epsR =  0.21; % based on Elas 
    minWinf = 10^(-7); maxWinf = 10^8;


elseif i == 7 % Squids
    A0 = exp([3.15 2.83 3.47]);
    b = [0]; 
    ratio = Winf/0.01; % 
    epsR =  0.03; % based on fish
    minWinf =  24.0709; maxWinf = 500000;
end


% 
    h0 = A0(j)/(param.epsa*(param.f0 - param.fc));

% growth coefficient: -----------------------------------------------------
    A = param.epsa.*h0.*Winf.^(b(1)).*(f - param.fc);

% population growth rate: -------------------------------------------------
    r_max = A.*(1 - n) .* Winf.^(n - 1) .* ((1 - param.a) .* log(ratio) + log(epsR));
    
% Minimum resource level:  ------------------------------------------------
    R_star = param.fc * h0 * Winf.^(b(1) + n - param.q)...
        ./(param.gamma * (1 - param.fc) ); 
    
end