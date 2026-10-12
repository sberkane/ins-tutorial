function Xnew = se23_propagate_classical(X, omega, a, g, dt)
% Classical INS propagation (small-angle approximation)
% J_l = I, Q_l = 0.5 I

[R, v, p] = se23_extract(X);

% Rotation
Rnew = R * so3_exp(omega*dt);

% Velocity
vnew = v + g*dt + R * (a*dt);

% Position
pnew = p + v*dt + 0.5*g*dt^2 + R * (0.5*a*dt^2);

Xnew = Gamma(Rnew, vnew, pnew);
end
