function result = poem(param, result)

%
% Init
%
y0 = param.y0;
 
%
% Run:
%
[t,y] = ode23(@poem_deriv, 0:1:param.tEnd, y0, odeset('NonNegative', 1:length(y0)), param);
%
% Construct output:
%
result.y = y;
result.R = y(:,param.ixR); % resource
result.B = y(:,param.ixFish);
result.t = t;

result.Yield = result.B .* (ones(length(t),1)*param.F');

B = mean(result.B(end-40:end,:)); % result.B(end,:); % 
y = mean(y(end-40:end,:)); % y(end,:); % 

[result.f, result.mortpred, result.Eavail] = calcEncounter(y', param);
result.mort = result.mortpred(param.ixFish)' + param.mort0 + param.F;
[result.v , result.nu] = calcNu(result.Eavail, result.mort, param);



