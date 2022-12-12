function param = baseparameters(WinfSquid)
% It would be great to have this model in a more flexible way: where we
% could specify the parameter physiology/prey pred interaction... with
% default value based on fish

%% TO DO
% Check changes imposed by changing the maturation size. 

%%
clear param
if nargin == 0
    WinfSquid = 3.5*10^3;
end

% PLot parameters: 
param.SpId = {'SmallPel','Mesopelagic','LargePel', 'Demersal', 'Squid'};
param.LWidth = [1.5, 1.5, 3, 3, 3];
param.Color = [0.07, 0.62, 1 ; ... 
    0.85, 0.33, 0.10; ...
    0.00, 0.45, 0.74; ...
    0, 0, 0 ;...
    1.00,0.55,0.00]; 

param.tEnd =  300;
param.ixR = [1 2 3 4]; % 4 resources: Small - large zoo small - large benthos 
param.w(param.ixR) = [2e-06  0.001 0.5e-03 0.25]; % weight lower limit
param.wc(param.ixR) = [2e-06*sqrt(500) 0.001*sqrt(500) 0.5e-03*sqrt(250000) 0.25*sqrt(500)]; % weight central size
param.wu(param.ixR) = [0.001 0.5  125  125]; % weight upper limit (Small benthos cover large benthos) 

%
% stages: -----------------------------------------------------------------
%   
param.nstage = 6; % number of stages predator use 3, 6, 9, etc (prey = 2/3)
param.nsize  = param.nstage + 1; % 
param.sfish = 0.001; % smallest size fish (all fish)
param.lsfish = 250; % Winf for small forage and meso
param.lfish = 125000; % largest size fish (only predator)
param.smat = 250*0.28; %   0.5; %   weight at maturity forage/meso
param.lmat = 125000*0.28; %   250; % weight at maturity predators
% cephalopods: 
param.sceph = 0.01; % smallest size all cephalopod
param.lceph = WinfSquid; % 3.5e+03; % largest size allcephalopod 5.5902e+03
param.SmallCeph = 100; % small cephalopod speceis
param.sizes = logspace(log10(param.sfish), log10(param.lfish), param.nsize); % All the size on log scale 
param.sizes_Ceph = logspace(log10(param.sceph), log10(param.lceph), param.nsize);
[~, param.maxsmall] = min(abs(param.sizes-250));

%
% Species: ----------------------------------------------------------------
% 
% Indices and weight classes small pelagics
param.w = [param.w param.sizes(1:param.maxsmall-1)];
param.ix1(1) = 4 + 1;
param.ix2(1) = 4 + (param.maxsmall-1);
% Indices and weight classes mesopelagics:
param.w = [param.w param.sizes(1:param.maxsmall-1)];
param.ix1(2) = param.ix2(1) + 1; 
param.ix2(2) = param.ix2(1) + (param.maxsmall-1);
% Indices and weight classes large pelagic:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(3) = param.ix2(2) + 1; 
param.ix2(3) = param.ix2(2) + (param.nsize-1);
% Indices and weight classes large demersal:
param.w = [param.w param.sizes(1:end-1)];
param.ix1(4) = param.ix2(3) + 1; 
param.ix2(4) = param.ix2(3) + (param.nsize-1);
% Indices and weight classes Squid:
param.w = [param.w param.sizes_Ceph(1:end-1)]; % intermediate size
param.ix1(5) = param.ix2(4) + 1; 
param.ix2(5) = param.ix2(4) + (param.nsize-1);

% Index for fish:  
param.nSpecies = length(param.ix1) ; %% add one group
param.ixFish = param.ix1(1):param.ix2(end); % 2 groups with 4 stages and 4 groups (incl. cephalopods) with 6 stages
param.NbSizeGrp = param.ixFish(end);

% Creating central and upper size limits from param.w (lower limit): with
% potential difference in size structure. 
param.stage = param.w(param.ix2)./param.w(param.ix2-1); % log distance between two size bin (for each group).
for i = 1:param.nSpecies 
    idx = param.ix1(i):param.ix2(i); 
    param.wc(idx) = param.w(idx).*sqrt(param.stage(i)); % central sizes
    param.wu(idx) = param.w(idx).*param.stage(i); % Upper sizes
end 

%
% predator prey preference ------------------------------------------------
%
param.beta(1:param.ixFish(end)) = 400; % optimal pred/prey mass ratio fish
param.beta(param.ix1(5):param.ix2(5)) = 50; % optimal pred/prey mass ratio Cephalopod 
param.sigma = 1.3;

% Size-preference matrix
thetaP = zeros(length(param.wc));
thetaP(5:end, :) = sqrt(pi/2)*param.sigma.* (...
 erf((log(param.wu) - log(param.wc(param.ixFish)./param.beta(param.ixFish))')./(sqrt(2)*param.sigma)) ... 
 - erf((log(param.w) - log(param.wc(param.ixFish)./param.beta(param.ixFish))')./(sqrt(2)*param.sigma)))...
 ./ log(param.wu./param.w);
param.prefer = thetaP; 

%
% Physiology: -------------------------------------------------------------
%
param.epsAssim = 0.7;                   % assimilation efficiency Fish
param.epst = 0.1;                       % efficiency of benthos community 
param.h = zeros(1, param.NbSizeGrp); param.h(1:end) = 20; % Value fish
param.h(param.ix1(5):param.ix2(5)) = 28/(param.epsAssim*(0.6 - 0.2)); % Value ceph
param.met = 0.2*param.h;  % maintenance costs, 20% of h
param.q = 0.8; % clearance rate exponent:
param.n = 3/4; % maximum consumption rate exponent
param.m = 0.825; % metabolic cost exponent 
param.gamma = 70; % factor for the max clearance rate (area per time)  
param.eRepro = repmat(0.01,param.nSpecies,1)';
param.A = param.epsAssim*param.h*(0.6 - 0.4); % here I assume a constant feeding level. 
param.mort0 = (0*param.ixFish'+.1);

%
% Temperature scaling on physiology
%
% Q10 =  1.88; % Assumption of similar Q10 for all the groups. 
param.Q10 = repelem(1.88, param.NbSizeGrp);
param.Q10m = repelem(1.88, param.NbSizeGrp);

% Physio: 
param.z = (param.w./param.wu)';
param.Cmax = (param.h.*param.wc.^param.n)./param.wc;
param.V = (param.gamma*param.wc.^param.q)./param.wc;
param.Mc = (param.met.*param.wc.^param.m)./param.wc;

% Fishing mortality:-------------------------------------------------------
% Traul Fishing: (small and large species)
param.Fi = [0.3, 0, 0.3, 0.3, 0.3]; % default fishing iphensity.
F = [];
% Fishing size preference: 
for i = 1:param.nSpecies
    F = [F, param.Fi(i) * (1+(param.wc(param.ix1(i):param.ix2(i))...
        ./(param.wu(param.ix2(i))*0.05)).^(-3)).^(-1)];
end
param.F = F'; 

% Maturation investment: --------------------------------------------------
[~,param.matstageS] = min(abs(param.sizes-param.smat));
[~,param.matstageL] = min(abs(param.sizes-param.lmat));

% Allometric function investment in growth: 
param.kappaS = 1-((1+(param.wc(param.ix1(1):param.ix2(1))./param.smat).^(-5)).^(-1).*...
    (param.wc(param.ix1(1):param.ix2(1))./param.lsfish).^(1-param.n));

param.kappaL = 1-((1+(param.wc(param.ix1(3):param.ix2(3))./param.lmat).^(-5)).^(-1).*...
    (param.wc(param.ix1(3):param.ix2(3))./param.lfish).^(1-param.n));
param.kappa = [param.kappaS param.kappaS param.kappaL param.kappaL [1 1 1 1 1 1]];

%
% Habitat: ----------------------------------------------------------------
%
param.photic = 150;
param.mesop = 250;
param.visual = 1.5; % scalar; >1 visual predation primarily during the day, = 1 equal day and night
param.S2P = 0.6; % predation from Squid to pelagics. 

%
% Vertical distribution: --------------------------------------------------
%
sigmad = 10; % width of initial distribution
tau = 10;  % increase in width
param.sigmap = sigmad + tau*log10(param.wc/param.wc(1)); % width for each size class

% first stages as juvenile/adult for predators
[~,param.ixjuv] = min(abs(param.sizes-param.smat));
[~,param.ixadult] = min(abs(param.sizes-param.lmat));

% For initial condition: 
param.B0 = 0*param.ixFish+.01;
% param.B0(param.ix1(5)-4:param.ix2(5)-4)=0; %put Squid to zero

end