function X = Gamma(R, v, p)
%GAMMA Embedding of (R,v,p) into SE_2(3) as a 5x5 matrix
%
% X = [ R  v  p
%       0  1  0
%       0  0  1 ]
%
% R: 3x3 rotation, v,p: 3x1 vectors

X = [R, v, p;
    0 0 0 1 0;
    0 0 0 0 1];
end
