function PlotEcosystem(depth, pprod, sqd)
    % Building app for Matlab app builder. The function first run the FEISTY
    % model, then produce a plot for ecosystem. 
    Ftsize = 12; 
    %
    %% Run FEISTY: 
    %
    
    % parameter for FEISTY 
    param = baseparameters();
    param = baseparam_depth(param, depth);
    param.K =  [pprod, pprod, 0, 0];   % g ww/m2
    param.y0 = [0.1*param.K 0.01*param.B0];  
    if depth < param.mesop
        param.y0(param.ix1(2):param.ix2(2))=0; % mesopelagics to zero
    end
    if sqd == 1 
        param.y0(param.ix1(5):param.ix2(5)) = 0;
    end
    
    % Run FEISTY: 
    ngroup = param.ix2(end);
    result = poem(param);
    % Average of the biomass:                                                     
    Bi = mean(result.y((end - 40):end,:)); % take values for the 40 last time steps
    
    % Calc detritus flux to Benthos
    [~, mortpred, ~] = calcEncounter(Bi, param); 
    F2B = (121 + 2.58*mortpred(1:2).*Bi(1:2))*param.martin*param.epst - param.r(3)'.*Bi(3) - mortpred(3)'.*Bi(3);

    %% Parameters for plot
    %
    % Calculate average depth day/night
    Av_depth = -(param.avlocDay + param.avlocNight)/2;
    % change a bit for visualisation: 
    Av_depth(param.ix1(1):param.ix2(1)) = Av_depth(param.ix1(1):param.ix2(1))+0.1*param.bottom;
    Av_depth(param.ix1(3):param.ix2(3)) = Av_depth(param.ix1(3):param.ix2(3))-0.1*param.bottom;
    if param.bottom < param.mesop
        Av_depth(param.ix1(5):param.ix2(5)) = Av_depth(param.ix1(5):param.ix2(5))+0.05*param.bottom;
        Av_depth(param.ix1(4):param.ix2(4)) = Av_depth(param.ix1(4):param.ix2(4))-0.05*param.bottom;
    end 
    
    % Marker size depends on biomass:--------------------------------------
    Msize = Bi;
    Msize(Msize == 0) = nan;
    Msize = log(Msize+1.01)*150;
    
    % Create flux from interaction: 
    % Creat coordinate for line between points: 
    coord_y = []; 
    coord_x = [];
    idx = [];
    for i = 1:ngroup
         idx(1,1:ngroup) = param.wc(1:ngroup);
         idx(2,1:ngroup) = param.wc(i);
         coord_x = [coord_x, idx];    
        
         idx(1,1:ngroup) = Av_depth(1:ngroup);
         idx(2,1:ngroup) = Av_depth(i);
         coord_y = [coord_y, idx];
    end 
    
    % Creat Linewidth: ----------------------------------------------------
    Theta = param.theta .* Bi .* Bi'; % flux equal the rate * the prey biomass (* 0 if pred = 0)
    Theta(3,1:2) = F2B;
    Theta = Theta(1:end); 
    Values = sort(Theta); 
    indx = find(Theta > Values(end-75)); % takes the x highest values of theta
    idxL = quantile(Theta(indx), [0.33, 0.66]); % Get quantiles of Theta distribution.  
   
    LineWdth = Theta; % give the appropriate size
    LineWdth(Theta >= idxL(2)) = 1.8;
    LineWdth(Theta < idxL(2) & LineWdth >= idxL(1)) = 0.9;
    LineWdth(Theta < idxL(1)) = 0.45;
    
    Linetype(Theta >= idxL(2)) = {'-'};
    Linetype(Theta < idxL(2) & LineWdth >= idxL(1)) = {'--'};
    Linetype(Theta < idxL(1)) = {':'};
    
    % Set up color:--------------------------------------------------------
    colorSet = zeros(ngroup, 3);
    colorSet(1:4,:) =  [0      0.5      0;
                 0.62   0.99     0.51; 
                 0.5    0.3      0;
                 0.5    0.3      0;];         
    
    for i = 1:param.nSpecies
        for j = param.ix1(i):param.ix2(i)
            colorSet(j, :) = param.Color(i,:);
        end
    end 
    
    Col = [];
    for i = 1:ngroup
        for j = 1:ngroup
        Col = [Col; colorSet(i,:)];
        end 
    end 
    
    % Drawing scatter plot: 
    p = plot(coord_x(:,indx) , coord_y(:,indx), 'k');
    
    for i = 1:length(indx)
        p(i).LineWidth = LineWdth(indx(i));
        p(i).Color = Col(indx(i), : );
        p(i).LineStyle = Linetype{indx(i)};
    end 
    
    hold on 
    scatter(param.wc, Av_depth, Msize, colorSet, 'filled')

    % Additional editting 
    set(gca,'BoxStyle','full','Color',...
    'none' ,'FontSize', Ftsize,'GridLineStyle',...
    'none','XMinorTick','on','XScale','log',... 
    'Ylim', [-1.10*param.bottom, max(Av_depth)+0.1*param.bottom], ...
    'Xlim', [10^(-5) 10^5], 'XTick', [10^(-3) 10^(1) 10^(5)], ...
    'YTick', [- param.bottom 0]);

end