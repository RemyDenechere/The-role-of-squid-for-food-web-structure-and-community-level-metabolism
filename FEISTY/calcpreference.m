function[theta, depthDay, depthNight, avlocDay, avlocNight] = calcpreference(param)
%
% Calculate overlap from depth distribution
%

% This part can be write in a more efficient way: matrix of Xloc for all
% the size class for day and night: Then just unsing it once in the main
% equation. 

xrange = linspace(0, param.bottom,(param.bottom+1)); % depths
% vertical migration depth: 
dvm = param.photic + 500; 
dvm(param.bottom < (param.photic + 500)) = param.bottom; % migration to bottom in intermediate habitats
dvm(param.bottom <= param.mesop) = 0; % no migration in shallow habitats

%
% zooplankton: ------------------------------------------------------------
%
% Night
xloc = 0; 
zp_n = (1./(sqrt(2*3.14*param.sigmap(param.ixR(1:2)).^2))) .* ...
        exp(-((xrange-xloc).^2' ./(2*param.sigmap(param.ixR(1:2)).^2)));
zp_n = zp_n*diag(1./sum(zp_n,1));

% Day (half at surface, half at dvm depth
xloc = dvm;
ix = param.ixR(1:2);zp_d = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc).^2' ./(2*param.sigmap(ix).^2)));
zp_d = zp_d*diag(1./sum(zp_d,1));

zp_d = (zp_n + zp_d)/2;

%
% benthos small and large: ------------------------------------------------
% (at bottom with width sigma)
xloc = param.bottom; 
bent = (1/(sqrt(2*3.14*10^2))) .* exp(-((xrange-xloc).^2/(2*10^2)));
bent = bent /sum(bent); 

bent_dn = [bent' bent']; 

%
% Small pelagic fish: -----------------------------------------------------
% (day + night)  
xloc = 0; 
ix = param.ix1(1):param.ix2(1);
spe = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc).^2' ./(2*param.sigmap(ix).^2)));

spel_dn = spe*diag(1./sum(spe,1));

%
% Meso pelagic: -----------------------------------------------------------
%
% night
mpel_n = spel_dn;

% day (at dvm)
xloc = dvm; 
ix = param.ix1(2):param.ix2(2);
mpe = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc).^2' ./(2*param.sigmap(ix).^2)));
mpel_d = mpe*diag(1./sum(mpe,1));

%
% large pelagic fish: -----------------------------------------------------
%
% night (surface)
xloc = 0; 
ix = param.ix1(3):param.ix2(3);
lpe = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc).^2' ./(2*param.sigmap(ix).^2)));
lpel_n = lpe*diag(1./sum(lpe,1));

% day (surface + dvm)
xloc = zeros(param.nstage,1)';
xloc(param.ixadult:end) = dvm;
ix = param.ix1(3):param.ix2(3);
lpe = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc').^2' ./(2*param.sigmap(ix).^2)));
lpel_d = lpe*diag(1./sum(lpe,1));
lpel_d = (lpel_d + lpel_n)/2; 

%
% demersal fish: ----------------------------------------------------------
%
% night
xloc = zeros(param.nstage,1)';
xloc(param.ixjuv:end) = param.bottom;
ix = param.ix1(4):param.ix2(4);
dem = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc').^2' ./(2*param.sigmap(ix).^2)));
dem_n = dem*diag(1./sum(dem,1));

% day
demmig = dvm;
demmig((param.bottom - dvm) >= 1200) = dvm + (param.bottom-dvm-1200);
demmig((param.bottom - dvm) >= 1500) = param.bottom;
xloc(param.ixadult:end) = demmig;

ix = param.ix1(4):param.ix2(4);
dem = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc').^2' ./(2*param.sigmap(ix).^2)));
dem_d = dem*diag(1./sum(dem,1));

% if shallower than euphotic depth, adult demersals feed across-habitats
if (param.bottom <= param.photic) 
   dem_d = (dem_d + dem_n)/2;
   dem_n = dem_d;
end

%
% Squids: -----------------------------------------------------------------
%
xloc = zeros(param.nstage,1)'; % Reset xloc = 0     
ix = param.ix1(5):param.ix2(5);
if param.bottom > param.mesop
   % AS Mesopelagic.  
   xloc(1:end) = dvm;
else  % In shallow water: AS demersal. 
   % large cephalopods 
   xloc = param.bottom;  
end
cph_d = (1./(sqrt(2*3.14*param.sigmap(ix).^2))) .* ...
        exp(-((xrange-xloc').^2' ./(2*param.sigmap(ix).^2)));
cph_d = cph_d*diag(1./sum(cph_d,1));
cph_n = lpel_n;

%
% calculate overlap during day and night: ---------------------------------
%
depthDay = [zp_d bent_dn spel_dn mpel_d lpel_d dem_d cph_d];
depthNight = [zp_n bent_dn spel_dn mpel_n lpel_n dem_n cph_n];
dayout = magic(param.ixFish(end)).*0;
nightout = magic(param.ixFish(end)).*0; 

for i = 1:param.ixFish(end)
   test = min(depthDay(1:end,i) , depthDay(1:end,1:end));
   dayout(1:end,i) = sum(test,1)';
   test =  min(depthNight(1:end,i), depthNight(1:end,1:end));
   nightout(1:end,i) = sum(test,1)'; 
end

% visual predation (pelagics) is good at light, bad in the dark) including
visualpred = [(param.ix1(1):param.ix2(1)) (param.ix1(3):param.ix2(3))]; % Spel Lpel 
dayout(visualpred, 1:(param.ix1(5))) = dayout(visualpred, 1:(param.ix1(5)))*param.visual;
nightout(visualpred, param.ix1(1):param.ix2(5)) = nightout(visualpred, param.ix1(1):param.ix2(5))*(2-param.visual); % Here we inplicitly assume that all prey have less pred from pel at night.
  
% Large pelagic predators limited vision in twilight zone during day
% Pelagics: 
pelpred = param.ix1(3):param.ix2(3); % Selects Lpel
pelpred = pelpred(param.ixadult:end); % Only adults
preytwi = param.ix1(2):param.ix2(2); % prey are Mesopelagics and bathypelagics
dayout(pelpred,preytwi) =  dayout(pelpred,preytwi)/param.visual*(2-param.visual);


% Squids: 
ceph = param.ix1(5):param.ix2(5);
if param.bottom > param.mesop % In deep regions Cepohaopods are hidding in the mesopelgaic regions at night (less predation)
                              % and are goin up at night (less predation freom visual
                              % predators). 
   dayout(visualpred, ceph) = dayout(visualpred, ceph)*0.75;
   nightout(visualpred, ceph) = nightout(visualpred, ceph)*0.75;
end 

% Average overlap during the whole day
location = (dayout+nightout).*0.5;

% calculate combined preference matrix   
theta = param.prefer.*location;

%
% change specific interactions -----------------------------------------
%

% benthivory is a strategy (larvae + pelagics do not eat benthos)
idx_be = [5:(param.ix1(4)+(param.ixjuv-2)), param.ix1(5):param.ix2(5)] ; % all larvae + pelagics (incl cepha). 
theta(idx_be, 3:4) = 0;   

% small demersals are less preyed on
idx_smd = param.ix1(4)+(param.ixjuv-1):param.ix1(4)+(param.ixadult-2);
theta(idx_be,idx_smd)= theta(idx_be,idx_smd)*.25;

% demersals do not eat much zooplankton
idx_dems = param.ix1(4)+(param.ixjuv-1):param.ix2(4);
theta(idx_dems,1:2) =  theta(idx_dems,1:2)*0;

% provide benefit to forage and mesopelagic fish (predator avoidance)
pred1 = param.ix1(3)+(param.ixadult-1):param.ix2(3); % large pelagic
pred2 = param.ix1(4)+(param.ixadult-1):param.ix2(4); % demersal

prey1 = param.ix1(1)+(param.ixjuv-1):param.ix2(1); % small pelagics
prey2 = param.ix1(2)+(param.ixjuv-1):param.ix2(2); % mesopel
idx_predat = [pred1 pred2];
idx_prey= [prey1 prey2];
theta(idx_predat,idx_prey) = theta(idx_predat,idx_prey)*0.5; % former value: 0.5


% Predation from Squid:
theta(ceph, (param.ix1(3)):(param.ix2(3))) = param.S2P* theta(ceph, (param.ix1(3)):(param.ix2(3)));
theta(ceph, (param.ix1(4)):(param.ix2(4)-4)) = param.S2P* theta(ceph, (param.ix1(4)):(param.ix2(4)-4));
theta(ceph, prey1) = param.S2P* theta(ceph, prey1);

% calculate center of vertical distribution during the day
[~,idi] = max(depthDay);
avlocDay = xrange(idi);

% calculate center of vertical distribution during the night
[~,idi] = max(depthNight);
avlocNight = xrange(idi);
 
  
