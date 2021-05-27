function OUT = avg(IN,varargin)
% OUT = avg(IN,dim)
%
% This function is identical to diff except that it calculates the
% average of adjacent elements instead of the difference. i.e., for
% dim=1, OUT = avg(IN,1) = (IN(2:end,:,:,:)+IN(1:(end-1),:,:,:))/2;
%
if (nargin<2)
    sz = size(IN);
    if (sz(1)>1)
        dim = 1;
    elseif (sz(2)>1)
        dim = 2;
    elseif (sz(3)>1)
        dim = 3;
    elseif (sz(4)>1)
        dim = 4;
    end
else
    dim = varargin{1};
end     
        
if (dim==1)
    OUT = (IN(2:end,:,:,:)+IN(1:(end-1),:,:,:))/2;
elseif (dim == 2)
    OUT = (IN(:,2:end,:,:)+IN(:,1:(end-1),:,:))/2;
elseif (dim == 3)
    OUT = (IN(:,:,2:end,:)+IN(:,:,1:(end-1),:))/2;
elseif (dim == 4)
    OUT = (IN(:,:,:,2:end)+IN(:,:,:,1:(end-1)))/2;
end

end
