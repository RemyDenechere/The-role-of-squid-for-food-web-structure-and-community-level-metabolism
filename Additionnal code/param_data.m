function [param] = param_data()
% Gathering all the data for Squid somatic growth rate A, egg size M_0, and
% A and M_0 data from Denechere et al 2022 AmNat for other groups 

%% Data for Squids: 
load("A_Squid.mat")
load("Egg_size_Squid.mat")
load("length_weight_Squid.mat")
load("lengthagerelationship.mat")

%% Length to weight coefficient: 
param.L2W.l = length_weight_Squid.Lengthcm;
param.L2W.w = length_weight_Squid.Weghtg;

%% weight at age: 
Tab = lengthagerelationship; 
Tab = Tab(~(isnan(Tab.Aged) | isnan(Tab.Weghtg) | Tab.Aged <= 0 | Tab.Weghtg == 0),:); % Delet rows containing Na or 0 for age and weight 

param.W2A.w = Tab.Weghtg;
param.W2A.a = Tab.Aged;
param.W2A.Winf = Tab.AsymptoticWeightg;
param.W2A.Sp = Tab.Species;
param.W2A.SpID = categories(Tab.Species) ; % list of species
param.W2A.lth = length(param.W2A.SpID);

%% A & Egg Squid: 
param.A.S = A_Squid.A;
param.AWinf.S = A_Squid.Winf;
param.S.S = A_Squid.Species;

param.w0.S = Egg_size_Squid.Hatching_weight_g;
param.Winf.S = Egg_size_Squid.Max_weight_g; 
param.Fit_s = exp(4.617)*param.Winf.S.^(0.9833);
%% Data from Denechere et al 2022:
load('Egg_size_data_Benthos.mat')
load('Egg_size_data_Elasmobranch.mat')
load('Egg_size_data_Teleost.mat')
load('Egg_size_data_Copepods.mat')
load('Egg_size_data_Mammal.mat')
load('Growth_data_Copepod.mat')
load('Growth_data_Bivalve.mat')
load('Growth_data_Elasmobranch.mat') 
load('Growth_data_Teleost.mat')

% Data from Denechere et al 2022

% conversion coefficient: 
param.B_c = 3.2*10^(-2); % Benthos length to weight coefficient g/cm^3 
param.B_eta = 16*10^(-2); % g, Benthos ratio between asymptotic weight (M) and weight at maturation
param.F_c = 0.01; % conversion coefficient from Length to weight g/cm^3 (see in Andersen 2019, book)
                  % F_c is used for both teleost and elasmobranches. 
param.F_eta = 0.28; % g, ratio between asymptotic size and size at maturation

% Individual growth rate A : ----------------------------------------------
% Bivalves:
Linf_B = Growth_data_Bivalve.Linf; % Asymptotic length bivalves (cm)
param.Winf_B = param.B_c*Linf_B.^(3); % conversion length (cm) to wet weigth (g)
K_B = Growth_data_Bivalve.K; % von Bertalanffy growth coef K 
% Conversion from K to A 
param.A_B = 3*K_B.*Linf_B.^(3/4)*param.B_c^(1/4)*param.B_eta^(-1/12); % g^(-1/4) year^(-1) with metabolic exponent n = 3/4 see Andersen 2019 (book) 
% fit of A with Winf for bivalves: 
param.FitAB = exp(-0.01941)*param.Winf_B.^(0.2149);


% Teleost parameters:
Linf_F = Growth_data_Teleost.Linf; % Linf = Asymptotic length (cm)
param.Winf_F = param.F_c*Linf_F.^(3); % Asymptotic weight (g) / conversion length to Wet weight
param.A_F = Growth_data_Teleost.A; % Somatic growth teleost g^(-1/4) year^(-1)
% of A with Winf for teleosts:
param.FitAF = exp(1.292)*param.Winf_F.^(0.04331); 


% Elasmobranch parameters: 
Linf_E = Growth_data_Elasmobranch.Linf; % Asymptotic length (cm)
param.Winf_E = param.F_c*Linf_E.^(3); % Asymptotic weight (g) / conversion length to Wet weight
K_E = Growth_data_Elasmobranch.K;  % von Bertalanffy growth coef K 
% Conversion from K to A Elas: 
param.A_E = 3*K_E.*Linf_E.^(3/4)*param.F_c^(1/4)*param.F_eta^(-1/12); % somatic growth: g^(-1/4) year^(-1)
% fit A with Winf:
param.FitAE = exp(1.706)*param.Winf_E.^(0.02404);


% Copepods 
param.ixAct = Growth_data_Copepod.Feeding == 'Active feeders' | Growth_data_Copepod.Feeding == 'Mixed feeders'; % selecte active and mixed feeders. 
param.ixPas = param.ixAct == 0;

param.A_C = Growth_data_Copepod.A; 
param.Winf_C = Growth_data_Copepod.Winf;
% fit model: 
param.FitAC_Act = param.Winf_C(~isnan(param.Winf_C(param.ixAct))).^(0)*6.8346; 
param.FitAC_Pass = param.Winf_C(~isnan(param.Winf_C(param.ixPas))).^(0)* exp(0.3652);


% Egg size strategy -------------------------------------------------------
% Get Linf [cm]: 
Linf_B = Egg_size_data_Benthos.L_i;
Linf_E = Egg_size_data_Elasmobranch.Linf;
Linf_T = Egg_size_data_Teleost.Linf;

% Get Winf: [g] conversion from length to weight
param.Winf_B2 = param.B_c*Linf_B.^(3); % 
param.Winf_E2 = param.F_c*Linf_E.^(3); % 
param.Winf_T2 = param.F_c*Linf_T.^(3); % 
param.Winf_C2 = Egg_size_data_Copepods.AdultSize*10^(-3); % mg to g
param.Winf_M2 = Egg_size_data_Mammal.AdultSize; % 

% get w0: [g]
param.w0_B = Egg_size_data_Benthos.Ww_0;
param.w0_E = Egg_size_data_Elasmobranch.w0;
param.w0_T = Egg_size_data_Teleost.w0;
param.w0_C = Egg_size_data_Copepods.ProgenySize*10^(-3); % mg to g
param.w0_M = Egg_size_data_Mammal.ProgenySize; %

% Fit: (power law) 
param.Fit_B = exp(14.28)*param.Winf_B2.^(1.07); 
param.Fit_E = exp(8.718)*param.Winf_E2.^(-0.2819); 
param.Fit_T = exp(7.30)*param.Winf_T2.^(0.9085);
param.Fit_C = exp(8.871)*param.Winf_C2.^(0.3256);  
param.Fit_M = exp(1.29)*param.Winf_M2.^(0.1294); 

%% param Rmax: 
param.a = 0.42; % physiological mortality year^(-1)
param.epsa = 0.6; % food assimilation efficiency for fish 
param.fc = 0.2; % critical feeding level
param.f0 = 0.6; % expected average feeding level for fish
param.q = 0.8; % Exponent for clearance rate
param.gamma = 1.9753e+03; % [g^(-q)L/yr]

%% Epsilon_R from elasmobranches
    
% r_max parameters for elasmobranches: 
A = 5.5; % somatic growth coeff
b = exp(1.681); % scaling exponent for somatic growth rate
ratio = 362.5060; % Adult:offspring size ratio
param.n = 3/4;
param.epsR_elasm = exp(b/(A*(1-param.n)) - (1-param.a)*log(ratio)) ;

end

