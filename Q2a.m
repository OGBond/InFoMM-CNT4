function [decVars, revenue,pricesOverTime] = Q2a(A,capVec,tInit,plotBool)

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

prLevels = 0:10:2000;
% Price levels (to be used as discretisations for each journey j)
K = length(prLevels);

pricesOffered = zeros(n,K);
for j = 1:n
    pricesOffered(j,:) = prLevels;
end
pricesOffered = pricesOffered';

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
%             if t >= 7
%                 dR = @(p,j) exp(-0.8*prEl(j)*(p/prices(j) - 1));
%                 else
%                 dR = @(p,j) exp(-prEl(j)*(p/prices(j) - 1));
%             end
            for k = 1:K
                objConst(k,j,i) = d(t)*dR(pricesOffered(k,j),j)*pricesOffered(k,j);
                ineqConstTemp(k,j,i) = d(t)*dR(pricesOffered(k,j),j)*A(l,j);
            end
        end
    end
    ineqConst(:,l) = ineqConstTemp(:);
end

for r = 1:tRange*n
    eqConst(r,K*r-(K-1):K*r) = 1;
end

[decVars,revenue] = linprog(-objConst(:)', ineqConst',capVec,eqConst,ones(tRange*n,1),...
    zeros(tRange*n*K,1),ones(tRange*n*K,1));

% Errors that were made:
% (1) Mistakenly assumed that because of the negative sign of the objective
% function (due to the problem actually being maximisation and linprog
% solving minimisation problems) that the bounds should be made negative
% and swapped
% (2) Originally made the price levels too small and there were not enough
% of them

% Deal with indices of each variable so that we can deduce each of the 
% time, route and price levels are used in the solution

decIndices = find(decVars==1);
possIndices_k = repmat(1:K,1,tRange*n);
possIndices_j = repmat(imresize(1:n, [1 n*K], 'nearest'),1,tRange);
possIndices_i = imresize(1:tRange, [1 tRange*n*K], 'nearest');

possIndices = [possIndices_k;possIndices_j;possIndices_i];

decVarIndices = possIndices(:,decIndices);

realQuantities = decVarIndices;

legs = [12, 13, 16, 23, 24, 32, 34, 35, 46, 52, 53]';
journeys = [12, 13, 124, 135, 16, 23, 24, 235, 246, 32, 34, 35, 346, 46, 52, 53, 534, 5346]';

for i = 1:length(decVarIndices)
    realQuantities(1,i) = prLevels(decVarIndices(1,i));
    realQuantities(2,i) = journeys(decVarIndices(2,i));
    realQuantities(3,i) = decVarIndices(3,i) + tInit - 1;
end
    
timeSeries = zeros(n,tRange);

for i = 1:tRange
    t = i + tInit - 1;
    for j = 1:n
        info = realQuantities(:,(realQuantities(3,:) == t) & (realQuantities(2,:)==journeys(j)));
        if length(info) == 3
           timeSeries(j,i) = info(1);
        end
    end
end
    
pricesOverTime = timeSeries;
   
revenue = -revenue;

legsIndicesInJourneys = [1,2,5,6,7,10,11,12,14,15,16];

demandOverTime = zeros(m,tRange);
for i = 1:tRange
    t = i - 1 + tInit;
    for l = 1:m
        demandOverTime(l,i) = d(t)*dR(pricesOverTime(i),legsIndicesInJourneys(l));
    end
end

cumDemandOverTime = cumsum(demandOverTime')';

if plotBool == 1
    subplot(1,2,1)
    plot(tInit:tFin,pricesOverTime')
    xlabel('Time after bookings open (days)')
    ylabel('Price of journey (£)')
    legend(cellstr(num2str(journeys)),'Location','eastoutside')
    subplot(1,2,2)
    plot(tInit:tFin,cumDemandOverTime')
    xlabel('Time after bookings open (days)')
    ylabel('Cumulative expected demand')
    legend(cellstr(num2str(legs)),'Location','eastoutside')
end

end