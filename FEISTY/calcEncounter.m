function [f, mortpred, Eavail, Enc] = calcEncounter(y, param)
%
% Consumption
%
Enc = param.V' .* param.theta.*y';

Encspecies = sum(Enc');
f = Encspecies./(param.Cmax+Encspecies); % correction for total prey consumption 
f(isnan(f)) = 0;
Eavail = param.Cmax.*param.epsAssim.*f - param.Mc;

%
% Mortality: from prey on i
mortpr =  ((param.Cmax .* param.V)' .* param.theta./ (param.Cmax + Encspecies)') .* y;
mortpred = sum(mortpr); % sum by column 
