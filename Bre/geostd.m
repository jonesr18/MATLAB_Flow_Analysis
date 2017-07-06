function out = geostd(M)

out = [];
for c = 1:size(M,2)
    V = M(:,c);
    out = [out exp(sqrt(sum(log(V./geomean(V)).^2)/length(V)))];
end

end