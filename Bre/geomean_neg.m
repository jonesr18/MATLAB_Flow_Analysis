function out = geomean_neg(data)

N = length(data);

pos_data = data(data>0);

neg_data = data(data<0);

if ~isempty(neg_data) && ~isempty(pos_data)
    out = geomean(pos_data)*length(pos_data)/N + geomean(-1.*neg_data)*length(neg_data)/N;
elseif ~isempty(neg_data)
    out = geomean(-1.*neg_data)*length(neg_data)/N;
elseif ~isempty(pos_data)
    out = geomean(pos_data)*length(pos_data)/N;
else
    out = NaN;
end
    


end