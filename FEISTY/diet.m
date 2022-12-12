function diet(depth, pprod,  squid)
% run model for diet plot 
% Squid = presence (1) or absence 0 of squid group

% Run model shallow + Squid  ----------------------------------------------
param = baseparameters();
param = baseparam_depth(param, depth); % param and depth (m)
param.K =  [pprod, pprod, 0, 0];  % Set up the secondary production (g ww/m2 )
param.y0 = [0.1*param.K 0.01*param.B0];
% Presence of Meso and midwater species depends on depth: 
if(param.bottom <= param.mesop) %
  param.y0(param.ix1(2):param.ix2(2))=0; % mesopelagics to zero
end
result = poem(param); 

w = param.wc;
y = result.y;
Bin = floor(0.8*length(y));
yend = mean(y(Bin:end,:));
ystage = param.ixFish(end);
ysmall = param.nstage - param.nstage*2/3;

[f, ~, ~, ~] = calcEncounter(yend', param);

bom = param.theta(5:ystage,:) .* mean(y(Bin:end,:)); 
fbom = f(5:ystage)' ./ sum(bom,2);
output = bom .* fbom;

colspec = [1 2 3 3  repmat(4,param.nstage*2/3,1)' ...
           repmat(5,param.nstage*2/3,1)'         ...
           repmat(6,param.nstage,1)'             ...
           repmat(7,param.nstage,1)'             ...
           repmat(8,param.nstage,1)'];
       
colorSet =  [0      0.5      0;
             0.62 0.99 0.51;
             0.5    0.3      0;
             param.Color];

SpId = {'SmallPel','Mesopelagic','LargePel', 'Demersal', 'Squid'};

if depth < 250 && squid == 1 % shallow system with squid: 
    nbrplt = 3; % 3 plots 
    idxplt = [1, 5, 4]; % (small pelagic, squid and demersal)
    dp = 'Shallow ';
    sq = 'With Squid';

elseif depth < 250 && squid == 0 % shallow system without squid:
    nbrplt = 2; % 2 plots 
    idxplt = [1, 4]; % (small pelagic and demersal) 
    dp = 'Shallow ';
    sq = 'Without Squid';

elseif depth >= 250 && squid == 1 % shallow system with squid:
    nbrplt = 3; % 3 plots (small pelagic, squid and large pelagic) 
    idxplt = [1, 5, 3];
    dp = 'Open Ocean ';
    sq = 'With Squid';

elseif depth >= 250 && squid == 0 % shallow system with squid:
    nbrplt = 2; % 3 plots (small pelagic and large pelagic) 
    idxplt = [1, 3];
    dp = 'Open Ocean ';
    sq = 'Without Squid';
end 

til = tiledlayout(1, nbrplt, Padding="compact");
for i = 1:nbrplt
    nexttile
    Sp_i = output(param.ix1(idxplt(i))-4:param.ix2(idxplt(i))-4,:);
    Sp_i = [Sp_i; zeros(ysmall, ystage)];
    H = bar(Sp_i, 'stacked');
     for j = 1:ystage
         H(j).FaceColor = colorSet(colspec(j),:);
         H(j).LineStyle = 'none';
     end
    ylim([0 1])
    title(SpId{idxplt(i)})
    ylabel('Feeding level') 
    set(gca,'XTick',[], 'FontSize', 14)
    xlabel('size-classes')
    Sp_i = [];
end
title(til, [dp, sq])
set(gcf, 'units', 'centimeters', 'position',[10, 10, 15, 5])

end