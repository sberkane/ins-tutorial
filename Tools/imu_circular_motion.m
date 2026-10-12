function [omega_b, a_b] = imu_circular_motion(t, V, r, g)
% Ideal IMU for planar circle in XY with yaw = (V/r) t.
% Returns:
%   omega_b: body angular rate (rad/s) in body frame
%   a_b    : specific force (m/s^2) in body frame

w=V/r;
yaw = w*t;
c = cos(yaw); s = sin(yaw);

R = [ c -s 0;
      s  c 0;
      0  0 1];

% World acceleration for p = [r sin(wt); r(1-cos(wt)); 0]
% => pddot = [-r w^2 sin; r w^2 cos; 0]
a_world = [ -r*w^2*sin(yaw);
             r*w^2*cos(yaw);
             0 ];

% Gyro: pure yaw rate about body z (same as world z here)
omega_b = [0;0;w];

% Accelerometer measures specific force: a_b = R' (a_world - g)
a_b = R' * (a_world - g);
end
