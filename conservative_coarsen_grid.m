function [Xout,Yout,Zout,Aout] = conservative_coarsen_grid(X,Y,Z,A,cfac);
% Coarsen grid of fluxes (Z) on lon-lat curvilinear grid (X,Y) with
% grid cell area (A) by factor cfac. X,Y,Z,A are 2D.
% e.g. Z is Wm-2, A is m2. cfac = 1 does nothing. cfac = 2 groups
% grid cells by 4, cfac = 3 groups grid cells by 9 etc.
% points are dropped at the end if the grid cell sizes don't match.


[xL,yL] = size(X);
xL = floor(xL/cfac)*cfac;
yL = floor(yL/cfac)*cfac;
xLout = xL/cfac;
yLout = yL/cfac;

X = X(1:xL,1:yL);
Y = Y(1:xL,1:yL);
A = A(1:xL,1:yL);
Z = Z(1:xL,1:yL).*A;
Z(isnan(Z)) = 0;

Xout = zeros(xLout,yLout);
Yout = Xout;
Zout = Xout;
Aout = Xout;

for xi=1:xLout
    if (mod(xi,50) == 0)
        ['Doing ' num2str(xi) ' of ' num2str(xLout)]
    end
    xinds = [0:(cfac-1)]+cfac*(xi-1)+1;
    for yi=1:yLout
        
        yinds = [0:(cfac-1)]+cfac*(yi-1)+1;
        
        Xout(xi,yi) = mean(mean(X(xinds,yinds)));
        Yout(xi,yi) = mean(mean(Y(xinds,yinds)));
        Aout(xi,yi) = sum(sum(A(xinds,yinds)));
        Zout(xi,yi) = sum(sum(Z(xinds,yinds)))/Aout(xi,yi);
    end
end

Zout(Zout==0) = NaN;

end
