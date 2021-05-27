function out = monmean(in,dim,weight)
% This function produces a weighted mean using the weight vector
% "weight" along the dimension dim.

    sz = size(in);
    if (length(sz)==2)
        sz = [sz(1) sz(2) 1 1];
    elseif (length(sz)==3)
        sz = [sz(1) sz(2) sz(3) 1];
    end
    
    total = sum(weight);

    if (dim == 1)
        out = sum(in.*repmat(weight(:),[1 sz(2:4)]),dim)/total;
    elseif (dim == 2)
        out = sum(in.*repmat(weight(:)',[sz(1) 1 sz(3:4)]),dim)/total;
    elseif (dim == 3)
        out = sum(in.*repmat(permute(weight(:),[3 2 1]),[sz(1:2) 1 sz(4)]),dim)/total;
    elseif (dim == 2)
        out = sum(in.*repmat(permute(weight(:),[4 3 2 1]),[sz(1:3) 1]),dim)/total;
    end
end

    
