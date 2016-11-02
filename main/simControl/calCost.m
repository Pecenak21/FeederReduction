function c = calCost(dat, costWeight)
% J = sum( Loss + costWeight * Power purchased)
c = sum(abs(dat.subLoss)) + costWeight * sum(real(dat.subPower));
end