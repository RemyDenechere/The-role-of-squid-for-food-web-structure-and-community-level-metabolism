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
        param.y0(param.ix1(5):param.ix2(5)) = 0; % calculation without squid!
        
        % Simulation with squid: ------------------------
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
                'Color', param.Color(sp,:), 'HandleVisibility',  HV2{dp}); 
        end
    end 
    % total biomass system: 
    plot(Zoo_prod, sum(Bi'), '--k', 'LineWidth', 1.5, 'HandleVisibility', HV{dp})  
   
    % Labels and legends: 
    if dp ==1
        ylim([0 40])
        ylabel('Biomass [g m^{-2}]')
    else 
        ylim([0 30])
        leg = legend([param.SpId{Sp_plot}, {'Total'}] , 'Location', 'best', ...
            'Color','none', 'EdgeColor', 'none');
    end
    title([titlelab(dp), '. ',num2str(depth(dp)), ' m'])
    set(gca, 'FontSize', 11)
    hold off

end % Depth

xlabel(t, 'Zooplankton productivity [g m^{-2} yr^{-1}]')