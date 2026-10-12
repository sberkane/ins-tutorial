function wx = so3_hat(w)
%SO3_HAT Hat operator for so(3)
% w: 3x1, wx: 3x3 skew-symmetric

wx = [   0   -w(3)  w(2);
    w(3)   0    -w(1);
    -w(2)  w(1)   0 ];
end
