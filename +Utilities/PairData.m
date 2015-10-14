function pareData = PairData(order1, order2)
% pair two vector, and resort the order2 return the sort result 'pareData'
% order2(paraData) = order1;
assert(numel(order1)==numel(order2))
pareData = zeros(size(order1));
[~,Ia] = sort(order1);
[~,Ib] = sort(order2);
pareData(Ia) = Ib;
pareData = int32(pareData);
% check result
assert(all((order2(pareData) - order1)< 10e-8));
end