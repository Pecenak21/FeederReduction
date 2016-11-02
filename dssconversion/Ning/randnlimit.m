function x = randnlimit(mu, std, minVal, maxVal, varagin)

assert(mu>=minVal && mu<=maxVal);
assert(std>0);

x = mu + std*randn(1,varagin);
outsideRange = x<minVal | x>maxVal;
i=1;
while sum(outsideRange)>0
   x(outsideRange) = mu + std*randn(sum(outsideRange),1);
   outsideRange = x<minVal | x>maxVal;
   i=i+1
end
end