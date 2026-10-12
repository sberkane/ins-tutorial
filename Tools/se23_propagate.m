function [Xnext] = se23_propagate(X, omega, a, g, dt)
%SE23_PROPAGATE Exact discrete propagation on SE_2(3):
%   X_{k+1} = G_dt * Phi_dt(X_k) * M_k
%
% Inputs:
%   X     : 5x5 SE_2(3) state (R,v,p)
%   omega : 3x1 measured angular rate (rad/s) (bias-free here)
%   a     : 3x1 measured specific force (m/s^2) in body frame (bias-free)
%   g     : 3x1 gravity vector in inertial frame (e.g., [0;0;9.81])
%   dt    : timestep (s)
%
% Outputs:
%   Xnext : propagated state

[Rk, vk, pk] = se23_extract(X);

% --- 1) Gravity increment
G = Gamma( eye(3), g*dt, 0.5*g*dt*dt );

% --- 2) Kinematics increment
Phi = Gamma( Rk, vk, pk + vk*dt );

% --- 3) IMU increment (piecewise-constant omega,a over [t,t+dt])
Rbar = so3_exp(omega * dt);
J    = so3_left_jacobian(omega * dt);
Q    = so3_left_Q(omega * dt);
vbar = J * (a*dt);
pbar = Q * (a*dt*dt);

M = Gamma(Rbar, vbar, pbar);

% --- Full update on the group
Xnext = G * Phi * M;
end
