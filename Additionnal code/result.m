%%
%              *Results, Discussion and Perspective* 
% 
addpath("Data")
addpath("C:/Users/rden/OneDrive - Danmarks Tekniske Universitet/Ph.D/Project_1/Deriving-pop-scaling-rules-from-individual-level-metabolism-and-life-history-trait/data/data in mat")

%% Result data. 
% DESCRIPTIVE TEXT
param = param_data();
col = my_color();
tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact")
nexttile(1)
hold on 
plot(param.Winf_B, param.A_B, 'o', 'Color', col.yell, 'MarkerSize', 5, 'MarkerFaceColor', [0.99 0.94 0.67]) % Benthos data A %
plot(param.Winf_F, param.A_F, '.', 'Color', [0, 0.26,0.99], 'MarkerSize', 7) % Teleost data A
plot(param.Winf_E, param.A_E, 'o', 'Color', col.redl, 'MarkerSize', 5) % Elasmobranch data A
plot(param.Winf_C(param.ixAct), param.A_C(param.ixAct), 'o', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82]) % Copepods data A active feeders 
plot(param.Winf_C(param.ixPas), param.A_C(param.ixPas), 'd', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82]) % Copepods data A Passive feeders  
plot(param.AWinf.S, param.A.S, 's', 'Color', col.ora, 'MarkerSize', 5, 'MarkerFaceColor', col.oral)
% fit: 
plot(param.Winf_B, param.FitAB, '-', 'Color', col.yel, 'Linewidth', 1.5) % Benthos fit A
plot(param.Winf_F, param.FitAF, '-', 'Color', col.blu, 'Linewidth', 1.5) % teleost fit A 
plot(param.Winf_E, param.FitAE, '-', 'Color', col.red, 'Linewidth', 1.5) % Elasmobranch fit A
plot(param.Winf_C(param.ixAct), param.FitAC_Act, '-', 'Color', [0.49,0.18,0.56], 'Linewidth', 1.5) % Copepods fit A
plot(param.Winf_C(param.ixPas), param.FitAC_Pass, '-', 'Color', [0.49,0.18,0.56], 'Linewidth', 1.5) % Copepods fit A
plot(param.AWinf.S, zeros(1, length(param.AWinf.S))+ mean(param.A.S), '-', 'Color', col.ora, 'Linewidth', 1.5)
hold off
ylabel('A (g^{1/4}/yr)')
%title('A')
set(gca, 'FontSize', 10, 'YScale', 'log', 'XScale', 'log', 'Box', 'on')                                                                                                                                                                                           
legend('Bivalve', 'Teleost', 'Elas.', 'Cop. A.F. ', 'Cop. A.P.', ...
    'Squid', 'Location', 'Southeast', 'EdgeColor', 'none', 'Color', 'none', ...
    'Position',[0.215667497012981 0.598265035115849 0.256613752987019 0.118755798217487], 'NumColumns', 2)

% Relation between egg size and Linf 
% 
nexttile(3)
% Data: 
hold on 
plot(param.Winf_B2, param.Winf_B2./param.w0_B, 'o', 'Color', col.yell, 'MarkerSize', 5, 'MarkerFaceColor', [0.99 0.94 0.67])
plot(param.Winf.S, param.Winf.S./param.w0.S, 'o', 'Color', col.ora, 'MarkerSize', 5, 'MarkerFaceColor', col.oral)
plot(param.Winf_T2, param.Winf_T2./param.w0_T, '.', 'Color', [0, 0.26,0.99], 'MarkerSize', 7)
plot(param.Winf_E2, param.Winf_E2./param.w0_E, 'o', 'Color', col.redl, 'MarkerSize', 5)
plot(param.Winf_C2, param.Winf_C2./param.w0_C, 'o', 'Color', [0.49,0.18,0.56], 'MarkerSize', 5, 'MarkerFaceColor', [0.78,0.6, 0.82])
plot(param.Winf_M2, param.Winf_M2./param.w0_M, 'o', 'Color', [0.51,0.66,0.31], 'MarkerSize', 5, 'MarkerFaceColor', [0.51,0.66,0.31]+ 0.3)
% Fit: 
plot(param.Winf_B2, param.Fit_B, '-', 'Color', col.yel, 'LineWidth', 1.5)
plot(param.Winf_T2, param.Fit_T, '-', 'Color', col.blu, 'LineWidth', 1.5)
plot(param.Winf_E2, param.Fit_E, '-', 'Color', col.red,  'LineWidth', 1.5)
plot(param.Winf_C2, param.Fit_C, '-', 'Color', [0.49,0.18,0.56],  'LineWidth', 1.5)
plot(param.Winf_M2, param.Fit_M, '-', 'Color', [0.51,0.66,0.31],  'LineWidth', 1.5)
plot(param.Winf.S, param.Fit_s, '-', 'Color', col.ora,  'LineWidth', 1.5)
hold off

xlabel('Asymptotic weight, M_{\infty} (g)')
ylabel('M_{\infty}/M_0 (g)')
%title('B')
legend('Bivalve', 'Teleost', 'Elas.', 'Cop.', 'Mammal', 'Location', 'best', ...
    'EdgeColor', 'none', 'Color', 'none', 'Box', 'on', ...
    'Position',[0.0622402355832951 0.308686044572019 0.169312166276748 0.18783068270595])
set(gca, 'FontSize', 10, 'YScale', 'log', 'XScale', 'log', 'Box', 'on')   

% Rmax simulation: 
color = [0.77,0.45,0.51; 0.00,0.45,0.74 ; 0.85,0.33,0.10 ; 0.93,0.69,0.13; 0.49,0.18,0.56 ; 0.49,0.18,0.56; 0.92,0.42,0.00]; % set the color 
color_light = [0.77,0.45,0.51; 0.41,0.76,0.99 ; 1.00,0.60,0.43 ; 1.00,0.82,0.39; 1.00,0.82,0.39]; 
% blue / red / Jaune / violet

[Winf, ~] = Grid(10^(-7), 10^8); % asymptotic size range g 
f = [0.3, 0.6, 1]'; % feeding level

nexttile(2, [2,1])
for i = 1:7
    [~ ,r_max, minWinf, maxWinf] = Pop_growth_rate(Winf, f(2), param, i, 1); % calculation rmax
    if i == 1 
        indx = [];
    else 
        indx = Winf >= minWinf & Winf <= maxWinf; 
    end

    % 90% confidence calculation:  
    [~ ,r_max_inf, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 2); % calculation conf. interval low boundary
    [~ ,r_max_sup, ~, ~] = Pop_growth_rate(Winf, f(2), param, i, 3); % calculation conf. interval upper boundary
    
    % Limitation for copepods and elasmobranchs
    if i == 3 % Elasmobranchs
        idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 1;
    elseif i == 5 | i == 6 % Copepods
        idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0 & Winf < 0.1;
    elseif i == 1 
        idx = Winf == Winf;
    else % Others
        idx = r_max_inf > 0 & r_max_sup > 0 & Winf > 0;
    end 

    if i == 6 
        plot(Winf(indx), r_max(indx), '',  Winf(idx), r_max(idx)', '--',...
        'Color', color(i,:), 'Linewidth', 1.5 )
    else 
    
    % Confident interval: based on variation of A from data collected:
    ciplot(r_max_inf(idx), r_max_sup(idx), Winf(idx), color(i,:), 0.25, 'none')
    hold on 
    plot(Winf(indx), r_max(indx), '-',  Winf(idx), r_max(idx)', '--',...
        'Color', color(i,:), 'Linewidth', 1.5 )
    end 
end

plot(Winf, 10^(0)*Winf.^(param.n-1), 'k--', 'LineWidth', 1.5);
ylabel('Maximum population growth rate r_{max} [yr^{-1}]')
xlabel('Adult size, M [g]')
%title('C')
set(gca, 'FontSize', 10, 'Xscale', 'log', 'Yscale', 'log', 'Xlim', [9.59252353156652e-08 101912836.96271], ...
    'Ylim', [0.00185088935383643 336.189436874531], 'Box', 'on')

save_pdf(gcf, 'Fig_Res1_Data_Rmax_Rstar', 20, 10)
