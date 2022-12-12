function dydt = poem_deriv(t, y, param)
y(y<0) = 0;
R = y(param.ixR);
B = y(param.ixFish);
%
% Feeding
%
[~, mortpred, Eavail] = calcEncounter(y, param);
%
% Fish:
%
ixFish = param.ixFish;

% Total mortality:
mort = mortpred(ixFish)' + param.mort0 + param.F ;         % param.vuln(ixFish)' .*                  

% Flux out of the size group:
v = (Eavail(param.ixFish))';  %./param.wc(param.ixFish))';
vplus = max(0,v);
gamma = (param.kappa'.*vplus - mort) ./ ...
    (1 - param.z(param.ixFish).^(1-mort./(param.kappa'.*vplus)) );
Fout = gamma.*B;
Repro = (1-param.kappa').*vplus.*B;                                        % A fraction 1-k of the available energy is
                                                                           % invested in reproduction.                                                                
% Flux into the size group(where Fout in last stage stays is included in repro):
for i = 1:param.nSpecies
    ix = (param.ix1(i):param.ix2(i)) - length(R);                          % select Size class in Species
    ixPrev = [ix(end) ix(1:(end-1))];           
    Fin(ix) = Fout(ixPrev);
    % Reproduction:
    Fin(ix(1)) = param.eRepro(i)*(Fin(ix(1)) + sum(Repro(ix)));                  
end
 
dBdt = Fin' - Fout + (v - mort).*B - Repro ;  

%
% Resource
%
dRdt(1:2) = param.r(1:2)'.*(param.K(1:2)'-R(1:2)) - mortpred(1:2)'.*R(1:2);

% detritus flux out: 
F_D_out = 121 + 2.58*(mortpred(param.ixR(1:2))*R(1:2));
dRdt(3) = F_D_out*param.martin*param.epst - param.r(3)'.*R(3) - mortpred(3)'.*R(3);
dRdt(4) = 0;

dydt = [dRdt'; dBdt];

