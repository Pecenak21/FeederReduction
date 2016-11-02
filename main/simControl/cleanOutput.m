function d = cleanOutput(d)
fn = fieldnames(d);
if d.k < d.niter
    k = d.k;
    for j = 1:length(fn)
        i = fn{j};
        if size(d.(i),1) < d.niter, continue; end;
        switch ndims(d.(i))
            case 2
                d.(i) = d.(i)(1:k,:);
            case 3
                d.(i) = d.(i)(1:k,:,:);
            otherwise
        end
    end
end
end