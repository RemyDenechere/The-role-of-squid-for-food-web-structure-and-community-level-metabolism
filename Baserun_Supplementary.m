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

%% To Do: 
% Sensitivity analysis: 
% - f, fc, in Rmax and Rstar calculation 
% - q and Gamme 
% - epsR, epsa : same effect 

%% set up
clear all 
close all   
addpath("Data\", "Data\Data from Denechere et al AmNat\", "FEISTY\", "Additionnal code\")

col = my_color();     % Colors for the different functional groups. 
param = param_data(); % Fish and Squid data and parameters for the rmax 
                      % and R* simulations. 
Save_Figures = true;
if Save_Figures
    mkdir Fig
end   

%% SupplementA Fig. SA1: length to weight coefficient: 
figure(1)

loglog(param.L2W.l, param.L2W.w, 'o', 'Color', 'k' , 'MarkerSize', 8, 'LineWidth', 1, ...
   'MarkerFaceColor', [0.75 0.75 0.75])
hold on 
plot(param.L2W.l, exp(-1.10666) *param.L2W.l .^ 2.21668, 'k', 'LineWidth', 1.5)
hold off 
xlabel('Squid length (cm)')
ylabel('Squid weight (g)')
set(gca, 'FontSize', 11)

%% SupplementA Fig. SA2: Weight at age cephalopod species:
figure(2)

param = param_data(); 
n = 3/4; % 2/3;
Winf = zeros(1, 11);

% Storage variable: 
Fit_para_nvar = zeros(2,param.W2A.lth);
Fit_para_cons = zeros(1,param.W2A.lth);

subplot(4,3,i)
t = tiledlayout(4, 3, "TileSpacing", "compact", "Padding", "compact");

for i = 1:param.W2A.lth

    idx = param.W2A.Sp == param.W2A.SpID{i};
    Winf(i) = max(param.W2A.Winf(idx));

    % Preparing data
    [x, y] = prepareCurveData(log(param.W2A.a(idx)), log(param.W2A.w(idx)) );
    % x = log(param.W2A.a(idx)); y = log(param.W2A.w(idx));
   
    % fitting n as variable: ----------------------------------------------
    ft = fittype( 'poly1' ); % Set up fittype
    [fitresult, ~] = fit(x, y, ft);

    % fitting n as const % ------------------------------------------------
    % Set up fittype and options:
    if n == 2/3
        ft_cons = fittype( 'a + x*3', 'independent', 'x', 'dependent', 'y' );
    elseif n == 3/4
        ft_cons = fittype( 'a + x*4', 'independent', 'x', 'dependent', 'y' );
    end 
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = 0.05; 
    [fitresult_cons, ~] = fit(x, y, ft_cons, opts); % Fit model to data.

        % get fitting parameters: 
    Fit_para_nvar(:,i) = [fitresult.p1 ; fitresult.p2]; 
    Fit_para_cons(i) = fitresult_cons.a; 

    % Plot data :----------------------------------------------------------
    nexttile
    loglog(param.W2A.a(idx), param.W2A.w(idx), 'o', 'color' , col.blul,...
        'MarkerSize', 2, 'MarkerFaceColor', col.blul);
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
xlabel(t, 'Age, t (days)')
ylabel(t, 'Body mass, m (g)'); 
legend('Data', 'Fit n var.', 'Fit n = 3/4', 'Location', 'northwest', 'Color', 'none')

% convert from [g^1-n/day] in per year: 
A_var = A_var * 365;
A_cons = A_cons * 365;

% Plot n for the species. 
% subplot(4,3,3*4)
% hist(n_var)
% ylabel('frequence'); xlabel('Metabolic exponent, $n$', 'Interpreter','latex') 

if n ==2/3
    disp(['Somatic growth rate, evaluated with n = 2/3: ' num2str(geomean(A_cons))])
elseif n == 3/4
    disp(['Somatic growth rate, evaluated with n = 3/4: ' num2str(geomean(A_cons))])
end

if Save_Figures
    save_graph(gcf, 'pdf', [ 'Fig/' 'Supplement_fit_w_age'], 16, 20)
end 

%% Calcul fit A and Winf: 

% fitting n_varying: ------------------------------------------------------
[xData, yData] = prepareCurveData(log(Winf), log(A_cons));
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft );

% assuming that A does not vary with winf use geommean of A 
% Set up fittype and options.
ft = fittype( 'a* x^0', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.636743444294958;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


% clear A_cons

%% SupplementA Fig. SA3: length at age with temperature:
figure
load('lengthagerelationship.mat') 
Tab = lengthagerelationship; 
Tab = Tab(~(isnan(Tab.Aged) | isnan(Tab.Weghtg) | isnan(Tab.TemperatureC) |...
    Tab.Aged <= 0 | Tab.Weghtg == 0),:); % Delet rows containing Na or 0 for age, weight and temp 
sp = categories(Tab.Species);            % list of species   
L = length(sp);                          % nbr of species                                    

t_layout = tiledlayout(4, 2, "TileSpacing", "compact");
t_layout.TileIndexing = "columnmajor";

k = 0;  % counter to 0
for i = 1:L % loop species
    indx_sp = Tab.Species == sp{i};
    Tab2 =  Tab(indx_sp ,:);                       % subtable per species    
    temp = unique(Tab2.TemperatureC);              % get the temperature experiments
    if isempty(temp)                               % next species if no temp
        continue
    end 
    ntemp = length(temp);                          % number of temperatures
    disp([sp{i}, ' ' ,num2str(temp'), '°C']);      % display species and temps

    for t = 1:ntemp  % loop temperature               
        k = k+1;                                   % counter + 1
        temp_val(k) = temp(t);                     % save all temperatures
        indx_temp = Tab2.TemperatureC == temp(t);  % Index tempr
        Tab3 =  Tab2(indx_temp ,:);                % subtable per species    
        
              
        % fitting n as variable: ----------------------------------------------
        x = log(Tab3.Aged); y = log(Tab3.Weghtg);
        [x, y] = prepareCurveData( x, y );                      % Preparing data  
        ft = fittype( 'poly1' );                                % Set up fittype
        [fitresult, ~] = fit(x, y, ft );                        % get fitting parameters
        Fit_para_nvar(:,k) = [fitresult.p1 ; fitresult.p2];     % save fitting param    

        % fitting n as Constant: ----------------------------------------------
        ft_cons = fittype( 'a + x*3', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = 0.05; 
        [fitresult_cons, ~] = fit(x, y, ft_cons, opts); % Fit model to data.
        Fit_para_cons(k) = fitresult_cons.a;

        % Calc fit :-----------------------------------------------------------
        %! variable
        n_var(k) = (Fit_para_nvar(1,k)-1)./Fit_para_nvar(1,k);  % Calc metabolic expoenent n
        b = (1-n_var(1, k)); % b = (1-n)
        A_var(k) = (exp(Fit_para_nvar(2,k)).^b)./b;             % Calc Somatic growth A
        a = linspace(min(Tab3.Aged), max(Tab3.Aged), 100);      % vector age
        y_var = ((1-n_var(k)).*A_var(k).*a).^(1./(1-n_var(k))); % Calc weight at age

        %! Constant: 
        n = 2/3;
       
        A_cons(k) = (exp(Fit_para_cons(k)).^(1-n))./(1-n); % asuming that n is constant
        y_cons = ((1-n)*A_cons(k)*a).^(1/(1-n)); % n const. 
        
        % Plot:
        nexttile
        hold on 
        plot(Tab3.Aged, Tab3.Weghtg, '.')
        plot(a, y_var, 'r', a, y_cons, '--b')
        hold off
        title({sp{i}, [num2str(temp(t)), ' °C']})      
    end 
end
% title(t_layout,'Somatic growth rate per species and temp')
xlabel(t_layout,'Age (d)')
ylabel(t_layout,'Weight (g WW)')

%% SupplementA Fig. SA4: Somatic growth rate with temperature:
figure
A_cons_peryr = A_cons * 365;
[xData, yData] = prepareCurveData(temp_val, A_cons_peryr);
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

figure(4) % A Vs Temp
plot( fitresult, xData, yData )
legend('Data', 'Fit', 'Location','best')
set(gca, 'XScale', 'linear')
xlabel('Temperature °C')
ylabel('Somatic growth rate (g^{1-n} yr^{-1})')


%% Supplement SC1: Ratio Cephalopods vs pelagic in ecopath models. 
load("Data\Dataecopathsceph.mat")
subplot(1,2,1)
plot(Dataecopathsceph.Fulx_Zpp_B_per_t, Dataecopathsceph.Ratio, 'o', 'Color', 'k' ,...
    'Linewidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', [0.75 0.75 0.75])
xlabel('Zooplankton productivity (g m-2 yr-1))')
ylabel('Ratio ceph vs. pelagics')
set(gca, 'FontSize', 11)
title('A')

subplot(1,2,2)
plot(Dataecopathsceph.Depth, Dataecopathsceph.Ratio, 'o', 'Color', 'k' ,...
    'Linewidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', [0.75 0.75 0.75])
xlabel('Depth (m)')
ylabel('Ratio ceph vs. pelagics (%gm^{-2})')
set(gca, 'FontSize', 11)
title('B')


%% Supplement SD1: Effect of squid efficiency to prey pelagic fish

% Set up range of parameters: ---------------------------------------------
figure
depth = [50 2000]; %                                     
Zoo_prod = 100;  % Vector for Zoo productivity 
param = baseparameters();
Spred = linspace(0, 1, 30);
Bi = zeros(length(Spred), param.nSpecies);
param.K =  [Zoo_prod, Zoo_prod, 0, 0];   % g ww/m2

tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact")

for dp = 1:length(depth)
       
    for i = 1:length(Spred)
        param.S2P = Spred(i);                       % Reset predation from squid on fish
        param = baseparam_depth(param, depth(dp));  % Recalcul predation matrix

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
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
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
            plot(Spred(2:end), Bi(2:end, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:)); 
        end
    end
    plot(Spred(2:end), sum(Bi(2:end, :)'), '--k', 'LineWidth', 1.5 )
    % Replot Mesopelagic for better visualisation: 
    sp = 2;
    plot(Spred(2:end), Bi(2:end, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:), 'HandleVisibility', 'off')
    hold off
    xlabel('Predation by Squid on Fish')
    if dp == 1
        ylabel('Biomass (g m-2) ')
    end 
    title([num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
end
legend([param.SpId{Sp_plot}, {'Tot'}] , 'Location', 'best', 'Color','none', ...
    'EdgeColor', 'none')


% -------------------------------------------------------------------------
% Effect of Fish predation intensity on squid (Squid as prey,
% top-down effect of fish on squid) 
% -------------------------------------------------------------------------

% Set up range of parameters: ---------------------------------------------
depth = [50 2000]; %                                     
Zoo_prod = 100;  % Vector for Zoo productivity 
param = baseparameters();
Fpred = linspace(0, 1, 30);
Bi = zeros(length(Spred), param.nSpecies);
param.K =  [Zoo_prod, Zoo_prod, 0, 0];   % g ww/m2

for dp = 1:length(depth)
       
    for i = 1:length(Spred)
        param.F2S = Fpred(i);                       % Reset predation from Fish on Squid
        param = baseparam_depth(param, depth(dp));  % Recalcul predation matrix

        if i == 1
            param.y0 = [0.1*param.K 0.01*param.B0];           
        else % start from previous final state. 
            param.y0 = results.y(end, :) + [0.1*param.K param.B0];
        end
        
        if(param.bottom <= param.mesop) %
            param.y0(param.ix1(2):param.ix2(2))=0; % mesopelagics to zero
        end
        
        results = poem(param);                                                    
        yend = results.y((end - 20):end,:); % take values for the 40 last time steps
        % sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
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
            plot(Spred(2:end), Bi(2:end, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:)); 
        end
    end
    plot(Spred(2:end), sum(Bi(2:end, :)'), '--k', 'LineWidth', 1.5 )
    sp = 2;
    plot(Spred(2:end), Bi(2:end, sp)', 'LineWidth', param.LWidth(sp), ... 
                'Color', param.Color(sp,:), 'HandleVisibility', 'off')
    hold off
    xlabel('Predation by Pelagic on Squid')
    if dp == 1
        ylabel('Biomass (g m-2) ')
    end 
    title([num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
end
%
if Save_Figures
    save_graph(gcf, 'pdf', [ 'Fig/' 'Supp_SensAnalSquidPred'], 16, 14)
end 

%% Supplement SE1: Turnover rate of energy of mature organism

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
    set(gca, 'FontSize', 11, 'YScale', 'log', 'YLim', [10^(-3) 25], ...
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
if Save_Figures
    save_graph(gcf, 'pdf', [ 'Fig/' 'Supp_Flux2SSB'], 16, 8)
end 

%%  Sensitivity analysis zooplankton productivity for varying parameters
% We resimulated fig 4 A and B, with changing paramters for squid to
% analyse the effect of unknown Squid parameters on the main results of
% this study: We simulated 2 depths for +- 10% of the value of the paramter used 
% Parameters simulated: 
% - f, fc, in Rmax and Rst

%% Supplement SG1: Sensitivity analysis Sigma
%                    Zooplankton productivity
%                    Shallow 

figure
% set up range of parameters: ---------------------------------------------
depth = 50;  % 
Zoo_prod = linspace(5, 150, 10); % vector for Zoo productivity 
param = baseparameters();
titlelab = ['A', 'B']; 
Sigma = [1, 400];

t = tiledlayout(3, 2, TileSpacing='compact', Padding= 'compact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Sigma 

for sig = 1:length(Sigma)
    param = baseparameters(3.5*10^3, Sigma(sig));
    dp = 1;
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
        yend = results.y((0.8*param.tEnd):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((0.8*param.tEnd):end, :);
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
                'Color', param.Color(sp,:)); 
        end
    end 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5)   % total biomass system with squid!
    
    % Labels and legends: 
    if dp ==1
        % ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(sig), '. \beta = ',num2str(Sigma(sig))])
    set(gca, 'FontSize', 11)
    hold off
end % Depth


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               A 
titlelab = ['C', 'D'];
A = [18  28];
for a_i = 1:length(A)
    param = baseparameters(3.5*10^3, 50, A(a_i));
    dp = 1;
    Bi = zeros(length(Zoo_prod), param.nSpecies);
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
        yend = results.y((0.8*param.tEnd):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((0.8*param.tEnd):end, :);
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
                'Color', param.Color(sp,:)); 
        end
    end 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5)   % total biomass system with squid!
    
    % Labels and legends: 
    if dp ==1
        % ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(a_i), '. A = ',num2str(A(a_i)), 'g^{1-n} yr^{-1}'])
    set(gca, 'FontSize', 11)
    hold off
end % A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Winf 
titlelab = ['E', 'F'];
Winf = [250  125000];
for wf = 1:length(Winf)
    param = baseparameters(Winf(wf));
    dp = 1;
    Bi = zeros(length(Zoo_prod), param.nSpecies);
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
        yend = results.y((0.8*param.tEnd):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((0.8*param.tEnd):end, :);
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
                'Color', param.Color(sp,:)); 
        end
    end 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5)   % total biomass system with squid!
    
    % Labels and legends: 
    if dp ==1
        % ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(wf), '. M_\infty = ', num2str(Winf(wf)), 'g'])
    set(gca, 'FontSize', 11)
    hold off
end % A
title(t, ['Shelf system, Depth = ', num2str(depth), ' m'])
xlabel(t, 'Zooplankton productivity [g m^{-2} yr^{-1}]')

if Save_Figures
    save_graph(gcf, 'pdf', 'Fig/Supp_SensAnalisis_shelfsystem', 16, 20)
end 

%% Supplement SG2: Sensitivity analysis Sigma
%                    Zooplankton productivity
%                    Deep ocean 

figure
% set up range of parameters: ---------------------------------------------
depth = 2000;  % m
Zoo_prod = linspace(5, 150, 10); % vector for Zoo productivity 
param = baseparameters();
titlelab = ['A', 'B']; 
Sigma = [1, 400];

t = tiledlayout(3, 2, TileSpacing='compact', Padding= 'compact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Sigma 

for sig = 1:length(Sigma)
    param = baseparameters(3.5*10^3, Sigma(sig));
    dp = 1;
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
        yend = results.y((0.8*param.tEnd):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((0.8*param.tEnd):end, :);
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
                'Color', param.Color(sp,:)); 
        end
    end 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5)   % total biomass system with squid!
    
    % Labels and legends: 
    if dp ==1
        % ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(sig), '. \beta = ',num2str(Sigma(sig))])
    set(gca, 'FontSize', 11)
    hold off
end % Depth


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               A 
titlelab = ['C', 'D'];
A = [18  28];
for a_i = 1:length(A)
    param = baseparameters(3.5*10^3, 50, A(a_i));
    dp = 1;
    Bi = zeros(length(Zoo_prod), param.nSpecies);
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
        yend = results.y((0.8*param.tEnd):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((0.8*param.tEnd):end, :);
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
                'Color', param.Color(sp,:)); 
        end
    end 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5)   % total biomass system with squid!
    
    % Labels and legends: 
    if dp ==1
        % ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(a_i), '. A = ',num2str(A(a_i)), 'g^{1-n} yr^{-1}'])
    set(gca, 'FontSize', 11)
    hold off
end % A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Winf 
titlelab = ['E', 'F'];
Winf = [250  125000];
for wf = 1:length(Winf)
    param = baseparameters(Winf(wf));
    dp = 1;
    Bi = zeros(length(Zoo_prod), param.nSpecies);
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
        yend = results.y((0.8*param.tEnd):end,:); % take values for the 40 last time steps

        % Sum of biomass per species average per time per size class per species :
        for ii = 1:param.nSpecies 
            Bi(i, ii) =  sum(mean(yend(:,param.ix1(ii):param.ix2(ii))));
        end 

        % Simulation without squid: --------------------
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        results = poem(param); 
        yend = results.y((0.8*param.tEnd):end, :);
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
                'Color', param.Color(sp,:)); 
        end
    end 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5)   % total biomass system with squid!
    
    % Labels and legends: 
    if dp ==1
        % ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(wf), '. M_\infty = ', num2str(Winf(wf)), 'g'])
    set(gca, 'FontSize', 11)
    hold off
end % A
xlabel(t, 'Zooplankton productivity [g m^{-2} yr^{-1}]')
title(t, ['Open Ocean, Depth = ', num2str(depth), ' m'])

if Save_Figures
    save_graph(gcf, 'pdf', 'Fig/Supp_SensAnalisis_OpenOcean', 16, 20)
end 