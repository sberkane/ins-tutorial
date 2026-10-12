function X = se23_exp(xi)
%SE23_EXP Exponential map for SE_2(3)
% xi = [omega; v; p] in R^9
%
% X = exp(xi^) = Gamma(R, V, P)

omega = xi(1:3);
v     = xi(4:6);
p     = xi(7:9);

% Rotation
R = so3_exp(omega);

% Left Jacobian of SO(3)
J = so3_left_jacobian(omega);

% Velocity and position
V = J * v;
P = J * p;

X = Gamma(R, V, P);
end
