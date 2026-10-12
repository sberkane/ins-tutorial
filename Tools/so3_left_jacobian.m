function J = so3_left_jacobian(w)
%SO3_LEFT_JACOBIAN Left Jacobian J_l(w) for SO(3)
% w: 3x1 rotation vector

theta = norm(w);
W = so3_hat(w);

if theta < 1e-3
    % Evaluate the coefficients without subtracting nearly equal numbers.
    theta2 = theta*theta;
    A = 1/2 - theta2/24 + theta2*theta2/720;
    B = 1/6 - theta2/120 + theta2*theta2/5040;
else
    A = (1 - cos(theta))/theta^2;
    B = (theta - sin(theta))/theta^3;
end
J = eye(3) + A*W + B*(W*W);
end
