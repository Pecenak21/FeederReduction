% function of disaggregating load
function IndiLoad = DisaggregLoad(AggreLoad,NumIndiLoad,std) % add paras RatingIndiLoad and BusNum
            for i = 1:1:length(AggreLoad)
%                 x = AggreLoad(i)/NumIndiLoad + std * randn(1,NumIndiLoad-1);
                x = randnlimit(AggreLoad(i), std, 0, 1, NumIndiLoad);
                X = x./sum(x);
                IndiLoad(i,:) = AggreLoad(i)*X;
%                 for j = 1:1:NumIndiLoad
%                     IndiLoad(i,j) = x(j);
%                 end
            end
end