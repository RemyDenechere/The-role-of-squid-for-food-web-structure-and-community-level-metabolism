function [Fout, Fin] = Get_Fout(y, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
       % set up variables:
        y(y<0) = 0;
        R = y(param.ixR)';
        B = y(param.ixFish)';

        % Calcul of flux out: 
        ixFish = param.ixFish;
        [~, mortpred, Eavail] = calcEncounter(y, param);

        % Total mortality:
        mort = mortpred(ixFish)' + param.mort0 + param.F;                  
        
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
end