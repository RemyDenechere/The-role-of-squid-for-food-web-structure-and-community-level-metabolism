%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       Base run supplementary Squid                      %
%                             Rémy Denéchère                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Supplementary Materials 'The role of squid for food web structure 
% and community-level metabolism' RémY Denéchère. 
% Parameters and their units are described in the 'param_data' for the
% rmax and R* rate data, and in 'baseparameter' for the FEISTY framework. 

%% set up
clear all 
close all   
addpath("Data\", "Data\Data from Denechere et al AmNat\", "FEISTY\")

col = my_color();     % Colors for the different functional groups. 
param = param_data(); % Fish and Squid data and parameters for the rmax 
                      % and R* simulations. 

%% SupplementB Fig. SB1: length to weight coefficient: 
figure(1)

loglog(param.L2W.l, param.L2W.w, 'o', 'Color', 'k' , 'MarkerSize', 8, 'LineWidth', 1, ...
   'MarkerFaceColor', [0.75 0.75 0.75])
hold on 
plot(param.L2W.l, exp(-1.10666) *param.L2W.l .^ 2.21668, 'k', 'LineWidth', 1.5)
hold off 
xlabel('Squid length (cm)')
ylabel('Squid weight (g)')
set(gca, 'FontSize', 12)

%% SupplementB Fig. SB2: Weight at age cephalopod species:
figure(2)

param = param_data(); 
n = 2/3;
Winf = zeros(1, 11);

% Storage variable: 
Fit_para_nvar = zeros(2,param.W2A.lth);
Fit_para_cons = zeros(1,param.W2A.lth);

for i = 1:param.W2A.lth

    idx = param.W2A.Sp == param.W2A.SpID{i};
    Winf(i) = max(param.W2A.Winf(idx));

    subplot(4,3,i)

    % Preparing data
    [x, y] = prepareCurveData(log(param.W2A.a(idx)), log(param.W2A.w(idx)) );
%     x = log(param.W2A.a(idx)); y = log(param.W2A.w(idx));
   
    % fitting n as variable: ----------------------------------------------
    ft = fittype( 'poly1' ); % Set up fittype
    [fitresult, ~] = fit(x, y, ft);

    % fitting n as const % ------------------------------------------------
        % Set up fittype and options:
    ft_cons = fittype( 'a + x*3', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = 0.05; 
    [fitresult_cons, ~] = fit(x, y, ft_cons, opts); % Fit model to data.

        % get fitting parameters: 
    Fit_para_nvar(:,i) = [fitresult.p1 ; fitresult.p2]; 
    Fit_para_cons(i) = fitresult_cons.a; 

    % Plot data :----------------------------------------------------------
    subplot(4,3,i)
    loglog(param.W2A.a(idx), param.W2A.w(idx), '.', 'color' , col.blul, 'MarkerSize', 8);
    xlabel('Age (day)'); ylabel( 'Weight (g)'); 
    title(param.W2A.SpID{i})
    hold on 
    
    % Plot fit :-----------------------------------------------------------
    n_var(i) = (Fit_para_nvar(1,i)-1)/Fit_para_nvar(1,i);
    b = (1-n_var(1,i)); % b = (1-n)
    A_var(i) = (exp(Fit_para_nvar(2,i)).^b)./b; % asuming that n is varying
    A_cons(i) = (exp(Fit_para_cons(i)).^(1-n))./(1-n); % asuming that n is constant
    
    % prediction using parameter evaluated: 
    a = linspace(min(param.W2A.a(idx)), max(param.W2A.a(idx)), 100); % vector age
    y_cons = ((1-n)*A_cons(i)*a).^(1/(1-n)); % n const. 
    y_var = ((1-n_var(i)).*A_var(i).*a).^(1./(1-n_var(i))); % n variying with sp
    
    % Add prediction on plot: 
    plot(a, y_var, 'r', a, y_cons, '--b')

end
legend('Data', 'Fit n var.', 'Fit n = 2/3', 'Location', 'northwest', 'Color', 'none')

% convert from [g^1-n/day] in per year: 
A_var = A_var * 365;
A_cons = A_cons * 365;

% Plot n for the species. 
subplot(4,3,3*4)
hist(n_var)
ylabel('frequence'); xlabel('Metabolic exponent, $n$', 'Interpreter','latex') 

%% Supplement SD1: Ratio Cephalopods vs pelagic in ecopath models. 
load("Data\Dataecopathsceph.mat")
subplot(1,2,1)
plot(Dataecopathsceph.Fulx_Zpp_B_per_t, Dataecopathsceph.Ratio, 'o', 'Color', 'k' ,...
    'Linewidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', [0.75 0.75 0.75])
xlabel('Zooplankton productivity (g m-2 yr-1))')
ylabel('Ratio ceph vs. pelagics')
set(gca, 'FontSize', 12)
title('A')

subplot(1,2,2)
plot(Dataecopathsceph.Depth, Dataecopathsceph.Ratio, 'o', 'Color', 'k' ,...
    'Linewidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', [0.75 0.75 0.75])
xlabel('Depth (m)')
ylabel('Ratio ceph vs. pelagics (%gm^{-2})')
set(gca, 'FontSize', 12)
title('B')


%% Supplement SE1: Effect of squid efficiency to prey pelagic fish

% Set up range of parameters: ---------------------------------------------
figure
depth = [50 2000]; %                                     
Zoo_prod = 100;  % Vector for Zoo productivity 
param = baseparameters();
Spred = linspace(0, 1, 30);
Bi = zeros(length(Spred), param.nSpecies);
param.K =  [Zoo_prod, Zoo_prod, 0, 0];   % g ww/m2

for dp = 1:length(depth)
       
    for i = 1:length(Spred)
        param = baseparam_depth(param, depth(dp));
        param.S2P = Spred(i);

        if i == 1
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) +[0.1*param.K param.B0];
        end
        
        if(param.bottom <= param.mesop) %
            param.y0(param.ix1(2):param.ix2(2))=0; % mesopelagics to zero
        end
        % param.y0(param.ix1(5):param.ix2(5)) = 0; 
        results = poem(param);                                                    
        yend = results.y((end - 20):end,:); % take values for the 40 last time steps
        % sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  mean(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 
    end 
    
    Sp_plot= [];
    % plot Biomass per species
    subplot(1,2,dp)
    hold on 
    for sp = 1:param.nSpecies
        if sum(Bi(:, sp)) == 0
        else % counting witch species are in the plot 
            Sp_plot = [Sp_plot, sp];                    
            idx = Bi(:, sp)>0; 
            plot(Spred(2:end), Bi(2:end, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:)); 
        end
    end
    plot(Spred(2:end), sum(Bi(2:end, :)'), '--k', 'LineWidth', 1.5 )
    hold off
    xlabel('Predation by Squid')
    if dp == 1
        ylabel('Biomass (g m-2) ')
    end 
    title([num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 12)
end
legend([param.SpId{Sp_plot}, {'Tot'}] , 'Location', 'northwest')

%% Supplement SF1: Turnover rate of energy of mature organism

figure
% set up range of parameters: ---------------------------------------------
depth = [50 2000];  % 250                                                
Zoo_prod = linspace(5, 150, 10); % vector for Zoo productivity 
param = baseparameters();
Fout = zeros(length(Zoo_prod), length(param.ixFish));
Fin = Fout; 
KAPPA = param.kappa; 
KAPPA(param.ix2 - 4) = 0;

for dp = 1:length(depth)
    param = baseparam_depth(param, depth(dp));

    for i = 1:length(Zoo_prod)
        param.K =  [Zoo_prod(i), Zoo_prod(i), 0, 0];   % g ww/m2

        if i == 1
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end
        
        if(param.bottom <= param.mesop) %
                param.y0(param.ix1(2):param.ix2(2)) = 0;
        end
        % param.y0(param.ix1(5):param.ix2(5))=0; 

        results = poem(param);                                                    
        y = mean(results.y((end - 40):end,:)); % take values for the 40 last time steps
        
        % Get Fin and Fout:
        [Fout(i,:), Fin(i,:)] = Get_Fout(y, param);
        Fout(i,:) = Fout(i,:).*y(5:end);
    end

    FOUT = zeros(length(Zoo_prod), param.nSpecies); 
    % Sum of the fluxes per group: 
    for ii = 1:param.nSpecies
        FOUT(:, ii) = sum((Fout(:, param.ix1(ii)-4:param.ix2(ii)-4).* ...
            [1-param.kappa(param.ix1(ii)-3:param.ix2(ii)-4), 1])'); 
        % multiply the flux out by the proportion of mature individual in
        % next size class. For the last size class Fout is going into
        % mature individual (1-kappa = 1). 
    end
    
    % plot Biomass per species
    subplot(1, length(depth), dp)
    set(gca, 'FontSize', 12, 'YScale', 'log', 'YLim', [10^(-3) 25], ...
        'XLim', [0 150])

    hold on 
    for sp = 1:param.nSpecies
        plot(Zoo_prod, FOUT(:, sp)', 'LineWidth', param.LWidth(sp), ... 
            'Color', param.Color(sp,:)); 

    end
    plot(Zoo_prod,  sum(FOUT'), '--', 'LineWidth', 3, 'Color', [0 0 0])
    hold off 
    if dp ==1
        ylabel({'Flux to SSB (g m-2 yr-1)'})
    else 
        legend([param.SpId{:}, {'Tot'}] , 'Location', 'northwest')
    end
    xlabel('Zooplankton productivity (g m-2 yr-1)')
    title([num2str(depth(dp)), 'm'])
end
%
% save_graph(gcf, 'pdf', 'Flux2SSB', 20, 10)

%% Supplement SF1: sensitivity analysis Winf

figure(3)
% set up range of parameters: ---------------------------------------------
depth = [50 2000];  % 250                                                
Zoo_prod = 40; % Zooplankton productivityg ww m-2 yr-1
WinfSquid = exp(linspace(log(250), log(125000), 10)); % Asymptotic size Squid ranging from 
                                    % small pelagic asymptotic size to large pelagic.
param = baseparameters(50);
titlelab = ['A', 'B']; 
tiledlayout(1, 2, TileSpacing='compact', Padding= 'compact')
Bi = zeros(length(WinfSquid), param.nSpecies);

for dp = 1:length(depth)
    for i = 1:length(WinfSquid)
        param = baseparameters(WinfSquid(i));
        param = baseparam_depth(param, depth(dp));
        param.K =  [Zoo_prod, Zoo_prod, 0, 0];

        if i == 1
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end
        
        % no mesopelagic in shallow system:
        if(param.bottom <= param.mesop) %
                param.y0(param.ix1(2):param.ix2(2)) = 0;
        end
        results = poem(param);                                                    
        yend = results.y((end - 40):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  mean(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 
    end 
    
    Sp_plot= [];
    % plot Biomass per species
    nexttile
    hold on 
    for sp = 1:param.nSpecies
        if sum(Bi(:, sp)) == 0
        else % counting witch species are in the plot 
            Sp_plot = [Sp_plot, sp];                    
            idx = Bi(:, sp)>0; 
            plot(WinfSquid, Bi(:, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:)); 
        end
    end
    plot(WinfSquid, sum(Bi'), '--k', 'LineWidth', 1.5 )
    hold off 
    if dp ==1
        ylabel('Biomass (g m-2)')
    else 
        legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none')
    end
    xlabel('Squid asymptotic size (g)')
    title([titlelab(dp), '. ',num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11, 'XScale', 'log')
end
% save_graph(gcf, 'pdf', 'SensAnaliWinf', 20, 10)
