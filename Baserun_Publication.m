%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       Base run publication Squid                        %
%                             Rémy Denéchère                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots for the manuscripts 'The role of squid for food web structure and
% community-level metabolism' RémY Denéchère. 
% Parameters and their units are described in the 'param_data' for the
% rmax and R* rate data, and in 'baseparameter' for the FEISTY framework. 

%%
clear all 
addpath("Data\", "Data\Data from Denechere et al AmNat\", "FEISTY\")

col = my_color();     % Colors for the different functional groups. 
param = param_data(); % Fish and Squid data and parameters for the rmax 
                      % and R* simulations. 

Save_Figures = true;
if Save_Figures
    mkdir Fig
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Figures publication                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fig 2: Day and night vertical distribution of fish and squid in FEISTY
figure(1)

p = baseparameters();
p = baseparam_depth(p, 500); % calcul vertical distribution for a bottom depth of 500 m

tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact')
TITLE = {'Shallow', 'Deep'};

for i = 1:5
    nexttile
    hold on
    plot(mean(p.depthDay(:,p.ix1(i):p.ix2(i))'), -(0:500), '-' , 'LineWidth', 2.5, ...
        'Color', p.Color(i, :) )
    plot(mean(p.depthNight(:,p.ix1(i):p.ix2(i))'), -(0:500), '--' , 'LineWidth', 2.5, ...
        'Color', p.Color(i, :))
    plot([min(mean(p.depthDay(:,p.ix1(i):p.ix2(i))')) ...
        max(mean(p.depthDay(:,p.ix1(i):p.ix2(i))'))], ...
        [-500, -500], 'k:')
    hold off
    if i == 1
        ylabel('Depth (m)')
    elseif i == 3
        xlabel('Average distribution')
        yticks([])
    elseif i == 4
        legend('Day','Night', 'Location','best')
        yticks([])
    else 
        yticks([])
    end
    set(gca, 'FontSize', 12)
    title(p.SpId{i}, 'FontSize', 12)
end 

if Save_Figures
    save_graph(gcf, 'pdf', [ 'Fig/' 'Vertical_pos_Squid'], 16, 8)
end 


%% Fig3: Data collected Squid and comparaison with Fish. 

figure(2)
tiledlayout(2, 2, TileSpacing="compact", Padding="compact")                % Set up the figure for plot

% Plot Growth parameter A: ------------------------------------------------
%
nexttile
hold on
plot(param.Winf_F, param.A_F, '.', 'Color', col.blul, 'MarkerSize', 7)     % Teleost data A vs. Winf (asymptotic size)
plot(param.AWinf.S, param.A.S, 's', 'Color', col.ora, 'MarkerSize', 7, ...
    'MarkerFaceColor', col.oral)                                           % Squid data A vs. Winf
plot(param.AWinf.S, 11.7635*param.AWinf.S.^0.114, '-', 'Color', ...
    col.ora,'Linewidth', 1.5)                                              % Teleost Fit data A vs. Winf
plot([min(param.Winf_F) max(param.Winf_F)],  [geomean(param.A_F(~isnan(param.A_F))), ...
    geomean(param.A_F(~isnan(param.A_F)))], ...
    'Color', col.blul,  'Linewidth', 1.5)                                  % Squid Fit data A vs. Winf
hold off
title('A')
ylabel('Somatic growth, A [g^{1-n}/yr]')
set(gca, 'FontSize', 12, 'XTick', [10^(0) 10^(2) 10^(4)], 'Xscale', 'log', ...
    'Yscale', 'log')                                                                                                                                                                                              

% Relation between offspring size and Winf: -------------------------------
% 
nexttile
% Data: 
hold on
loglog(param.Winf_T2, param.w0_T, '.', 'Color', col.blul, 'MarkerSize', 7) % Teleost data W0 vs. Winf 
plot(param.Winf.S,  param.w0.S, 's', 'Color', col.ora, 'MarkerSize', 7, ...
    'MarkerFaceColor', col.oral)                                           % Squid data W0 vs. Winf 
plot([min(param.Winf_T2) max(param.Winf_T2)], ...
    [geomean(param.w0_T(~isnan(param.w0_T))), geomean(param.w0_T(~isnan(param.w0_T)))], ...
    'Color', col.blul,  'Linewidth', 1.5)                                  % Teleost Fit data W0 vs. Winf 
plot([min(param.Winf.S) max(param.Winf.S)], [geomean(param.w0.S), geomean(param.w0.S)], ...
    'Color', col.ora,'Linewidth', 1.5)                                     % Squid Fit data W0 vs. Winf 
hold off
ylabel('Offspring size, M_0 [g]')
title('B')
set(gca, 'FontSize', 12, 'XTick', [10^(1) 10^(3) 10^(5)], 'Xscale', 'log', ...
    'Yscale', 'log')                                      


%
%% Rmax: -----------------------------------------------------------------
%
nexttile

[Winf, ~] = Grid(10^(-2), 10^7); % asymptotic size range g 
f = [0.3, 0.6, 1]'; % feeding level

hold on 
for i = [2, 7] % indixes for Fish and squid respectively 
    if i == 2
        Col = col.blu;
    else 
        Col = col.ora;
    end 
    [~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1); % calculation rmax
    % 90% confidence calculation:  
    [~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2);       % calculation conf. interval low boundary
    [~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3);       % calculation conf. interval upper boundary
    
    indx = Winf >= minWinf & Winf <= maxWinf; 
    idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0;

    % Confident interval: based on variation of A from data collected:
    ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), Col, 0.5, 'none')
    plot(Winf(indx), r_max(indx), '-',  Winf(idx), r_max(idx)', '--',...
        'Color', Col, 'Linewidth', 1.5 )
end 
hold off 
set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [10^(-1) 10^(7)], ...
    'YLim', [3*10^(-3) 2*10^(1)])
ylabel({'Max. pop. growth rate','r_{max} [yr^{-1}]'})
xlabel('Adult size, M [g]')
set(gca, 'FontSize', 12, 'XTick', [10^(-1) 10^(2) 10^(5)])
title('C')
box('on')

%
% R*:----------------------------------------------------------------------
%
nexttile

[R_sr_t , ~, ~, ~] = Pop_growth_rate(Winf, f(2), param, 2, 1); % calculation rmax
[R_sr_s , ~, ~, ~] = Pop_growth_rate(Winf, f(2), param, 7, 1); % calculation rmax
    
loglog(Winf, R_sr_t, 'Color', col.blu, 'LineWidth', 2)
hold on 
plot(Winf, R_sr_s, 'Color', col.ora , 'LineWidth', 2)
hold off
% Create doublearrow
annotation(gcf,'doublearrow',[0.689639005650495 0.689639005650495],...
    [0.422315705128205 0.259852300333368],'LineWidth',1);

% Create textbox
annotation(gcf,'textbox',...
    [0.694048675234286 0.281459001068376 0.157578928932379 0.107604879840938],...
    'String',{'Difference in','intraspecific','competition'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

xlim([0.0116039432421003 11603943.2421003]);
ylim([0.000714432566659917 0.00833606683109503]);

ylabel({'Minimum resource level', 'R^* [g]'})
xlabel('Adult size, M [g]')
set(gca, 'FontSize', 12, 'XTick', [10^(-1) 10^(2) 10^(5)])
legend('Teleost', 'Squid', 'Location', 'southwest')
title('D')

if Save_Figures
    save_graph(gcf, 'pdf', [ 'Fig/' 'Fig_Res1_Data_Rmax_Rstar'], 16, 13)
end 

%% Fig 4: Sensitivity analysis zooplankton productivity and 

%                    Zooplankton productivity

figure(3)
% set up range of parameters: ---------------------------------------------
depth = [50 2000];  % 250                                                
Zoo_prod = linspace(5, 150, 10); % vector for Zoo productivity 
param = baseparameters();
titlelab = ['A', 'B']; 
HV2= {'off', 'on'};
HV = {'on', 'off'};

tiledlayout(2, 2, TileSpacing='compact', Padding= 'compact')

for dp = 1:length(depth)
    Bi = zeros(length(Zoo_prod), param.nSpecies);
    Bi_nosquid = Bi;
    param = baseparam_depth(param, depth(dp));
    for i = 1:length(Zoo_prod)
        param.K =  [Zoo_prod(i), Zoo_prod(i), 0, 0];   % g ww/m2
        
        %! Initial condition:
        if i == 1 % Initial condition for 1st zooplankton conditions
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end

        % No mesopelagic in shallow waters! 
        if(param.bottom <= param.mesop) 
                param.y0(param.ix1(2):param.ix2(2)) = 0;
        end
        
        % Simulation with squid: ------------------------
        results = poem(param);                                                    
        yend = results.y((end - 40):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((end - 40): end, :);
        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi_nosquid(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 
    end
    
    % plot Biomass per species: -------------------------------------------
    Sp_plot= [];
    nexttile
    hold on 
    for sp = 1:param.nSpecies
        if sum(Bi(:, sp)) == 0
        else % Counting witch species are in the plot 
            Sp_plot = [Sp_plot, sp];                    
            idx = Bi(:, sp)>0; 
            plot(Zoo_prod, Bi(:, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:), 'HandleVisibility',  HV2{dp}); 
        end
    end 
    ciplot(sum(Bi'), sum(Bi_nosquid, 2), Zoo_prod, [0.75 0.75 0.75], 0.5, "none", "off")
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5, 'HandleVisibility', HV{dp})   % total biomass system with squid!
    plot(Zoo_prod, sum(Bi_nosquid, 2), '--', 'Color', [0.5 0.5 0.5] , ...
        'LineWidth', 1.5, 'HandleVisibility',  HV{dp});                             % total biomass system without squid!
    
    % Labels and legends: 
    if dp ==1
        ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
        legend('Tot. With Squid', 'Tot. Without Squid' , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(dp), '. ',num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
    hold off



end % Depth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Accumulation biomass 
%
% set up range of parameters: ---------------------------------------------
Met_cost = zeros(length(Zoo_prod), length(param.ixFish));
titlelab = ['C', 'D']; 

for dp = 1:length(depth)
    param = baseparam_depth(param, depth(dp));
    SUMF2mat = zeros(1, length(Zoo_prod));
    for i = 1:length(Zoo_prod)
        param.K =  [Zoo_prod(i), Zoo_prod(i), 0, 0];   % g ww/m2

        if i == 1
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end

        %! Simulation with squid: --------------------
        if(param.bottom <= param.mesop) %
                param.y0(param.ix1(2):param.ix2(2)) = 0;
        end
        results = poem(param);                                                    
        yend = mean(results.y((end - 40):end,:)); % take values for the 40 last time steps
        Met_cost(i,:) = param.Mc(param.ixFish) .* yend(param.ixFish);
    end
        
    % plot Biomass per species
    nexttile
    hold on 
    for sp = 1:param.nSpecies
        F2mat = sum(Met_cost(:, param.ix1(sp)-4:param.ix2(sp)-4)');
        if sum(F2mat) == 0
        else 
            ciplot(SUMF2mat, SUMF2mat + F2mat, Zoo_prod, param.Color(sp,:), 0.75, '-')
        end
        SUMF2mat = SUMF2mat + F2mat;
         
    end
    plot(Zoo_prod,  SUMF2mat, '--','LineWidth', 3, 'Color', [0 0 0])
    hold off 
    if dp ==1
        ylabel({'Metabolic loss, Mc_iB_i [g m^{-2} yr^{-1}]'})
    else 
        legend([param.SpId{:}, {'Total'}] , 'Location', 'northwest', ...
            'EdgeColor', 'none', 'Color', 'none')
    end
    xlabel({'Zooplankton productivity', '[g m^{-2} yr^{-1}]'})
    title([titlelab(dp), '. ', num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
end

if Save_Figures
    save_graph(gcf, 'pdf', [ 'Fig/' 'SensAnalZoop'], 16, 14)
end 

%% Fig 4bis: Sensitivity analysis zooplankton productivity without Squid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Zooplankton productivity

figure
% set up range of parameters: ---------------------------------------------
depth = [50 2000];  % 250                                                
Zoo_prod = linspace(5, 150, 10); % vector for Zoo productivity 
param = baseparameters();
titlelab = ['A', 'B']; 

tiledlayout(2, 2, TileSpacing='compact', Padding= 'compact')
Bi = zeros(length(Zoo_prod), param.nSpecies);

for dp = 1:length(depth)
    param = baseparam_depth(param, depth(dp));
    for i = 1:length(Zoo_prod)
        param.K =  [Zoo_prod(i), Zoo_prod(i), 0, 0];   % g ww/m2
        
        %! Initial condition:
        if i == 1 % Initial condition for 1st zooplankton conditions
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        % No mesopelagic in shallow waters! 
        if(param.bottom <= param.mesop) 
                param.y0(param.ix1(2):param.ix2(2)) = 0;
        end
        %param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        % Simulation:
        results = poem(param);                                                    
        yend = results.y((end - 40):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 
    end
    
    % plot Biomass per species: -------------------------------------------
    Sp_plot= [];
    nexttile
    hold on 
    for sp = 1:param.nSpecies
        if sum(Bi(:, sp)) == 0
        else % Counting witch species are in the plot 
            Sp_plot = [Sp_plot, sp];                    
            idx = Bi(:, sp)>0; 
            plot(Zoo_prod, Bi(:, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:)); 
        end
    end
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5 )  % total biomass system with squid!
    % plot(Zoo_prod, Bi, '-k' ,'LineWidth', 1.5);         % total biomass system without squid!
    
    % Labels and legends: 
    if dp ==1
        ylabel('Biomass [g m^{-2}]')
    else 
        legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none')
    end
    title([titlelab(dp), '. ',num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
    hold off 

end % Depth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Accumulation biomass 
%
% set up range of parameters: ---------------------------------------------
Met_cost = zeros(length(Zoo_prod), length(param.ixFish));
titlelab = ['C', 'D']; 

for dp = 1:length(depth)
    param = baseparam_depth(param, depth(dp));
    SUMF2mat = zeros(1, length(Zoo_prod));
    for i = 1:length(Zoo_prod)
        param.K =  [Zoo_prod(i), Zoo_prod(i), 0, 0];   % g ww/m2

        if i == 1
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        if(param.bottom <= param.mesop) %
                param.y0(param.ix1(2):param.ix2(2)) = 0;
        end
        results = poem(param);                                                    
        yend = mean(results.y((end - 40):end,:)); % take values for the 40 last time steps
        Met_cost(i,:) = param.Mc(param.ixFish) .* yend(param.ixFish);
    end
        
    % plot Biomass per species
    nexttile
    hold on 
    for sp = 1:param.nSpecies
        F2mat = sum(Met_cost(:, param.ix1(sp)-4:param.ix2(sp)-4)');
        if sum(F2mat) == 0
        else 
            ciplot(SUMF2mat, SUMF2mat + F2mat, Zoo_prod, param.Color(sp,:), 0.75, '-')
        end
        SUMF2mat = SUMF2mat + F2mat;
         
    end
    plot(Zoo_prod,  SUMF2mat, '--','LineWidth', 3, 'Color', [0 0 0])
    hold off 
    if dp ==1
        ylabel({'Metabolic loss, Mc_iB_i [g m^{-2} yr^{-1}]'})
    else 
        legend([param.SpId{:}, {'Total'}] , 'Location', 'northwest', ...
            'EdgeColor', 'none', 'Color', 'none')
    end
    xlabel({'Zooplankton productivity', '[g m^{-2} yr^{-1}]'})
    title([titlelab(dp), '. ', num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
end


%%  Plot ecosystem and feeding level with and without Squid: 
figure(4)

tiledlayout(2,2, Padding="compact", TileSpacing="compact");
Ftsize = 11;  % 14; 

nexttile
PlotEcosystem(50, 130, 0)
ylabel('Depth (m)')
title('A. With Squid')
set(gca, 'FontSize', Ftsize)

nexttile
PlotEcosystem(50, 130, 1)
title('B. Without Squid')
yticks([])
set(gca, 'FontSize', Ftsize)

nexttile
PlotEcosystem(2000, 130, 0)
% set(gca, 'YScale', 'log')
xlabel('weight (g)')
ylabel('Depth (m)')
title('C. With Squid')
set(gca, 'FontSize', Ftsize)
 
nexttile
PlotEcosystem(2000, 130, 1)
% set(gca, 'YScale', 'log')
xlabel('weight (g)')
yticks([])
title('D. Without Squid')
set(gca, 'FontSize', Ftsize)
%

%% plot diet ---------------------------------------------------------------
% Run in shallow (50m depth), productive (130 g/m3/yr), with squid:
figure(5)
diet(50, 130, 1)

% Plot Shallow (50m depth), productive (130 g/m3/yr), without  squid:
figure(6)
diet(50, 130, 0)

% Plot deep (2000m depth), productive (130 g/m3/yr), with squid:
figure(7)
diet(2000, 130, 1)

% Plot deep (2000m depth), productive (130 g/m3/yr), without squid:
figure(8)
diet(2000, 130, 0)








