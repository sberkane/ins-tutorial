function R = so3_exp(w)
%SO3_EXP Exponential map exp(so3_hat(w)) -> SO(3)
% w: 3x1 rotation vector (rad)

theta = norm(w);
W = so3_hat(w);

if theta < 1e-8
    % Series: exp(W) ≈ I + W + 1/2 W^2
    R = eye(3) + W + 0.5*(W*W);
else
    A = sin(theta)/theta;
    B = (1 - cos(theta))/theta^2;
    R = eye(3) + A*W + B*(W*W);
end
end
