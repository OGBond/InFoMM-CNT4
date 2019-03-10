function [optCtrl, revenue, lambda] = Q1a(A,capVec,tInit)

[m,n] = size(A);

% matA is the incidence matrix of the network, with legs on the row and
% journeys along the columns

% capVec is the vector of capacities on each of the legs, which must be a single
% vector with as many rows as A

% t is a time (positive integer) to be specified

tFin = 10;  % Total number of time steps

d = @(t) 300*0.5./(1 + exp(-0.5.*(t-5))); 
% Demand arrival rate as function of time

demandExp = zeros(n,1); 
% Storage space for expected demands in route-independent case

for t = tInit:tFin
    demandExp = demandExp + d(t);
end

prices = [220,220,400,250,200,230,200,200,200,200,230, ...
          120,150,150,200,150,160,230];
% Prices given in original table


[optCtrl,~,~,~,lambda] = ...
    linprog(-prices, [A;eye(n)], [capVec';demandExp],...
                     [],[], zeros(n,1));

revenue = prices*optCtrl;