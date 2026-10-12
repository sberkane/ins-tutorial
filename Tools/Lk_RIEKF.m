function Lk = Lk_RIEKF(X, ~, dt)
%LK_RIEKF First-order noise injection: dt times the continuous adjoint map.
% The gravity argument is retained for compatibility; it does not enter
% this first-order approximation.
[R, v, p] = se23_extract(X);
L11 = dt*R;   L12 = zeros(3);
L21 = dt*so3_hat(v)*R;   L22 = dt*R;
L31 = dt*so3_hat(p)*R;  L32 = zeros(3);

Lk = [L11 L12;
    L21 L22;
    L31 L32];
end
