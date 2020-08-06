function [ylin,b] = linfit(x,y)
xvec = [ones(length(x),1) x];
b = xvec\y;
ylin = xvec*b;
end
