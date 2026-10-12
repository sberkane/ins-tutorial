function [R, v, p] = circular_ground_truth(t, V, r)
%CIRCULAR_GROUND_TRUTH Analytic planar circular motion.
%
% Circle in XY plane, z = 0, yaw = w t.
% Body x-axis tangent to trajectory.
%
% p(t) = [ r sin(wt); r (1 - cos(wt)); 0 ]
% v(t) = [ V cos(wt); V sin(wt); 0 ]
% R(t) = Rz(wt)

w=V/r;
yaw = w*t;
c = cos(yaw); s = sin(yaw);

% Rotation: world <- body
R = [ c -s  0;
    s  c  0;
    0  0  1];

% Position
p = [ r*sin(yaw);
    r*(1 - cos(yaw));
    0 ];

% Velocity
v = [ V*cos(yaw);
    V*sin(yaw);
    0 ];
end
