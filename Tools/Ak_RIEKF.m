function Ak = Ak_RIEKF(g, dt)
%Discrete-time transition matrix for the RIEKF
Ak = [ ...
    eye(3),                zeros(3),         zeros(3);
    so3_hat(g)*dt,            eye(3),            zeros(3);
    0.5*so3_hat(g)*dt^2,      dt*eye(3),         eye(3) ];
end