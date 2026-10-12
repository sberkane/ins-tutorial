function Q = so3_left_Q(w)
%SO3_LEFT_Q Second-order left Jacobian Q_l(w) used for position increment:
% p_bar = Q_l(w) * a * dt^2, where w = omega*dt.
%
% Uses stable series for small angles.

theta = norm(w);
W = so3_hat(w);

if theta < 0.1
    % The fourth-order denominator makes the closed form sensitive to
    % cancellation. Series through theta^8 are accurate at this cutoff.
    theta2 = theta*theta;
    alpha = 1/6 + theta2*(-1/120 + theta2*(1/5040 + ...
        theta2*(-1/362880 + theta2/39916800)));
    beta = 1/24 + theta2*(-1/720 + theta2*(1/40320 + ...
        theta2*(-1/3628800 + theta2/479001600)));
else
    alpha = (theta - sin(theta))/theta^3;
    beta  = (theta^2 + 2*cos(theta) - 2) / (2*theta^4);
end
Q = 0.5*eye(3) + alpha*W + beta*W*W;
end
