% Calculates bid prices as per Q1b, with t + 1 = 2

BidPrices = zeros(1,18);

for j = 1:18;
    [o1, r1] = Q1a(A,c,2);
    [o2, r2] = Q1a(A,c-A(:,j),2);
    BidPrices(j) = r1 - r2;
end
