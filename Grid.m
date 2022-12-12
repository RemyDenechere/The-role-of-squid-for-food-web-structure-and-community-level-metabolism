function [w, dw] = Grid(W0, Winf)
%GRID creats a exponential grid for the size (from the egg size to the
% asymptotic size

w = W0; % egg size
j=1;
paramGrid = 0.01; % size step parameter

while w(j)< Winf  
    w(j+1) = w(j)*(1+paramGrid);
    j = j+1;
end
w(end) = Winf; % final value exactly the asymptotic size

% Size step values:
dw = w(2:end)- w(1:(end-1)); 
end

