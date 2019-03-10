function [decVars, revenue, lambda] = Q2a(A,capVec,tInit)

[m,n] = size(A);

% matA is the incidence matrix of the network, with legs on the row and
% journeys along the columns

% capVec is the vector of capacities on each of the legs, which must be a single
% vector with as many rows as A

% t is a time (positive integer) to be specified

tFin = 10;  % Total number of time steps

d = @(t) 300*0.5./(1 + exp(-0.5.*(t-5))); 
% Demand arrival rate as function of time

prEl = [1 1.2 2 1 0.8 1 0.9 2 1 1 2 2 2 1 1 2 1 1];
% Price elasticities for specific example

demandExp = zeros(tFin-tInit+1,n);
% Storage space for expected demands in route-independent case
    
prices = [220,220,400,250,200,230,200,200,200,200,230, ...
          120,150,150,200,150,160,230];
% Prices given in original table

demandResp = zeros(1,n);
% Storage space for expected demands in route-dependent case

dR = @(p,j) exp(-prEl(j)*(p/prices(j) - 1));
% Demand response to price as a function of p = p_{tj}, the price of
% journey j at time t

prLevels = [50,100,150,200,250,300,350,400];
% Price levels (to be used as discretisations for each journey j)
K = length(prLevels);

for j = 1:n
    prodLevel(j) = prLevels(min(find(prLevels>=prices(j))));
    % Approximates original prices by giving it the smallest price level
    % above it
end

tRange = tFin - tInit + 1;

objConst = zeros(K,n,tRange);
ineqConstTemp = zeros(K,n,tRange);
ineqConst = zeros(tRange*n*K,m);
eqConst = zeros(tRange*n,tRange*n*K);
% Coefficients of the decision variables z_{tjk} according to the
% expression on slide 42 of the lecture notes

for l = 1:m
    for j = 1:n
        for t = tInit:tFin
            i = t - tInit + 1;
            for k = 1:K
                objConst(k,j,i) = d(t)*dR(prodLevel(k),j)*prodLevel(k);
                ineqConstTemp(k,j,i) = d(t)*dR(prodLevel(k),j)*A(l,j);
            end
        end
    end
    ineqConst(:,l) = ineqConstTemp(:);
end

for r = 1:tRange*n
    eqConst(r,K*r-(K-1):K*r) = 1;
end

decVars = linprog(objConst(:)', ineqConst',capVec,eqConst,ones(tRange*n,1),...
    zeros(tRange*n*K,1),ones(tRange*n*K,1));
end