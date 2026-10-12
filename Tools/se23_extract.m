function [R,v,p] = se23_extract(X)
%SE23_EXTRACT Extract (R,v,p) from SE_2(3) matrix X (5x5)

R = X(1:3,1:3);
v = X(1:3,4);
p = X(1:3,5);
end
