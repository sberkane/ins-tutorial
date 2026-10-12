function result = loglinear_error(outputDir)
%LOGLINEAR_ERROR Compare exact invariant and local MEKF error models.
%
% Two nonlinear INS trajectories use the same noise-free IMU input and
% differ only in their initial condition. This deterministic experiment
% compares error dynamics without measurement updates.
%
%   Invariant coordinates:
%       xi_exact(t) = log( Xhat(t) X(t)^{-1} )^vee
%       xi_LL(t)    = exp(A_inv t) xi_0
%
%   MEKF coordinates used in the tutorial:
%       x_M,exact = [ log(Rhat R^T)^vee ; v-vhat ; p-phat ]
%       dot(x_M,lin) = A_M(t) x_M,lin,
%       A_M(t) = [ 0             0  0;
%                  (Rhat*a)^     0  0;
%                   0            I  0 ].
%
% The invariant model is exact in logarithmic coordinates for these
% group-affine dynamics. The MEKF model is only a first-order, local model;
% its discrepancy is therefore a model mismatch, not a navigation error.
%
% Required functions in Tools:
%   Gamma.m, se23_extract.m, se23_propagate.m,
%   imu_circular_motion.m, so3_exp.m, so3_hat.m
%

if nargin < 1, outputDir = []; end
[outputDir, cleanup] = experiment_setup(outputDir); %#ok<ASGLU>
figDir = outputDir;

%%-------------------------------------------
%  LaTeX / CSM plotting style

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FONT_AXES   = 22;
FONT_LABEL  = 25;
FONT_LEGEND = 22;
FONT_NOTE   = 18;
LINE_WIDTH  = 2;

bg_color   = [1.0, 0.9, 0.9];
axes_color = [1.0, 1.0, 1.0];
color_inv  = [0.00, 0.70, 1.00];  % blue: invariant/log-linear
color_mekf = [1.00, 0.50, 0.00];  % orange: MEKF

%%-------------------------------------------
%  Simulation settings

T  = 5.0;                         % duration [s]
dt = 0.01;                        % propagation step [s]
t  = 0:dt:T;

g = [0;0;9.81];                   % NED gravity [m/s^2]

% Both maneuvers start from the same physical state.
V0 = 5.0;                         % speed [m/s]
R0 = eye(3);
v0 = [V0;0;0];
p0 = zeros(3,1);
X0 = Gamma(R0,v0,p0);

r_circle = 12.5;                  % circle radius [m], omega_z = 0.4 rad/s

% Initial estimation error: attitude only in ordinary physical coordinates.
% Initially, vhat=v and phat=p.
theta0_deg = 20;
axis0 = [1;-1;0.5];
axis0 = axis0/norm(axis0);
Rerr0 = so3_exp(deg2rad(theta0_deg)*axis0);
Xhat0 = Gamma(Rerr0*R0,v0,p0);

% The right-invariant group error couples rotation and translation, so xi0
% generally has translational components even when vhat0=v0 and phat0=p0.
xi0 = se23_log_local(Xhat0*se23_inv_local(X0));

% Exact initial coordinates for the stated MEKF convention.
xM0 = mekf_error_coordinates(X0,Xhat0);

% Continuous-time exact log-linear matrix for the right-invariant error.
A_inv = [zeros(3),   zeros(3), zeros(3); ...
         so3_hat(g), zeros(3), zeros(3); ...
         zeros(3),   eye(3),   zeros(3)];

%%-------------------------------------------
%  Run the two maneuvers

resultStraight = run_case('straight',X0,Xhat0,xi0,xM0,A_inv, ...
                          t,dt,V0,r_circle,g);
resultCircle   = run_case('circular',X0,Xhat0,xi0,xM0,A_inv, ...
                          t,dt,V0,r_circle,g);

%%-------------------------------------------
%  Numerical checks / console summary

fprintf('\n===============================================================\n');
fprintf(' Exact log-linear error propagation example\n');
fprintf(' Initial attitude error: %.1f deg\n',theta0_deg);
fprintf(' Duration: %.2f s, dt = %.4f s\n',T,dt);
fprintf('===============================================================\n');

print_case_summary('Straight',resultStraight);
print_case_summary('Circular',resultCircle);

fprintf(['\nInterpretation:\n', ...
    '  * The invariant/log-linear discrepancy remains at numerical\n', ...
    '    precision for both trajectories.\n', ...
    '  * The MEKF discrepancy is a local-linearization model mismatch,\n', ...
    '    not a navigation-estimation error.\n\n']);

%%-------------------------------------------
%  Figure: quantitative model mismatch

fig = figure('Color',bg_color,'Position',[100 100 850 560]);

invPosMax = max([resultStraight.inv_pos_mismatch, ...
                 resultCircle.inv_pos_mismatch]);
invPosBound = one_significant_digit_upper_bound(invPosMax);

% ---- Position-coordinate mismatch ---------------------------------------
% Explicit margins prevent the axes box, legend, and terminal tick labels
% from crowding the exported page boundary. The normalized position leaves
% 9% top padding and 6% right padding.
ax = axes('Parent',fig,'Units','normalized', ...
          'Position',[0.11,0.15,0.83,0.76]);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');

hInv = plot(ax,t,resultStraight.inv_pos_mismatch,'-', ...
    'Color',color_inv,'LineWidth',LINE_WIDTH);
plot(ax,t,resultCircle.inv_pos_mismatch,'--', ...
    'Color',color_inv,'LineWidth',LINE_WIDTH);
hMekf = plot(ax,t,resultStraight.mekf_pos_mismatch,'-', ...
    'Color',color_mekf,'LineWidth',LINE_WIDTH);
plot(ax,t,resultCircle.mekf_pos_mismatch,'--', ...
    'Color',color_mekf,'LineWidth',LINE_WIDTH);

format_panel(ax,T,FONT_AXES,axes_color);
xlabel(ax,'$t$ (s)','FontSize',FONT_LABEL);
ylabel(ax,'Position mismatch (m)','FontSize',FONT_LABEL);
ylim(ax,[0,1.08*max([resultStraight.mekf_pos_mismatch, ...
                     resultCircle.mekf_pos_mismatch])]);

text(ax,0.97,0.08,latex_precision_note(invPosBound,'m'), ...
    'Units','normalized','FontSize',FONT_NOTE,'Color',color_inv, ...
    'BackgroundColor',axes_color,'Margin',3, ...
    'HorizontalAlignment','right','VerticalAlignment','bottom');

% Keep the primary legend focused on model/color.
lgd = legend(ax,[hInv,hMekf], ...
    {'Log-linear model','MEKF linearization'}, ...
    'FontSize',FONT_LEGEND,'Location','northwest', ...
    'Orientation','horizontal');
lgd.Box = 'off';

% MATLAB supports one legend per axes. Use a transparent overlay axes to
% provide a second legend with visible line-style samples.
axMotionKey = axes('Parent',fig,'Units','normalized', ...
    'Position',ax.Position,'Color','none','Visible','off', ...
    'HitTest','off','PickableParts','none');
hold(axMotionKey,'on');
hStraight = plot(axMotionKey,nan,nan,'-','Color',[0.15,0.15,0.15], ...
    'LineWidth',LINE_WIDTH);
hCircular = plot(axMotionKey,nan,nan,'--','Color',[0.15,0.15,0.15], ...
    'LineWidth',LINE_WIDTH);

lgdMotion = legend(axMotionKey,[hStraight,hCircular], ...
    {'Straight motion','Circular motion'}, ...
    'FontSize',FONT_LEGEND,'Orientation','horizontal');
lgdMotion.Box = 'off';
lgdMotion.Units = 'normalized';
lgdMotion.Position = [0.125,0.70,0.62,0.07];

% Vector export with the tutorial background preserved.
outFile = fullfile(figDir,'loglinear_vs_mekf.pdf');
drawnow;
exportgraphics(fig,outFile,'ContentType','vector', ...
               'BackgroundColor',bg_color,'Padding','figure');
fprintf('Figure exported to:\n  %s\n',outFile);
result = struct('straight',resultStraight,'circular',resultCircle,'t',t);
save(fullfile(outputDir,'loglinear_error.mat'),'result');
end

%% ========================================================================
%  Local functions

function out = run_case(caseName,X0,Xhat0,xi0,xM0,A_inv,t,dt,V,r,g)
%RUN_CASE Propagate the nonlinear states and both error-coordinate models.

N = numel(t);
X = X0;
Xhat = Xhat0;
xM_lin = xM0;

xi_exact_hist = zeros(9,N);
xi_ll_hist    = zeros(9,N);
xM_exact_hist = zeros(9,N);
xM_lin_hist   = zeros(9,N);

for k = 1:N
    % Both exact coordinate errors are formed from this same state pair.
    xi_exact = se23_log_local(Xhat*se23_inv_local(X));
    xM_exact = mekf_error_coordinates(X,Xhat);

    % Closed-form exact invariant prediction from the common initial error.
    xi_ll = expm(A_inv*t(k))*xi0;

    xi_exact_hist(:,k) = xi_exact;
    xi_ll_hist(:,k)    = xi_ll;
    xM_exact_hist(:,k) = xM_exact;
    xM_lin_hist(:,k)   = xM_lin;

    if k == N
        break;
    end

    % Bias-free IMU input held constant over [t_k,t_{k+1}].
    switch lower(caseName)
        case 'straight'
            % With NED gravity, R=I, omega=0, and a=-g imply dot(v)=0.
            omega = zeros(3,1);
            a = -g;
        case 'circular'
            [omega,a] = imu_circular_motion(t(k),V,r,g);
        otherwise
            error('Unknown motion case: %s',caseName);
    end

    % Evaluate the trajectory-dependent MEKF Jacobian at the interval
    % midpoint. The same zero-order-held IMU sample used by the nonlinear
    % propagation is used here. Midpoint evaluation plus expm(A_M*dt) is a
    % second-order accurate discretization and keeps its integration error
    % negligible relative to the local-linearization mismatch of interest.
    Xhat_mid = se23_propagate(Xhat,omega,a,g,0.5*dt);
    [Rhat_mid,~,~] = se23_extract(Xhat_mid);
    A_M = [zeros(3),            zeros(3), zeros(3); ...
           so3_hat(Rhat_mid*a), zeros(3), zeros(3); ...
           zeros(3),            eye(3),   zeros(3)];
    xM_lin = expm(A_M*dt)*xM_lin;

    % Both nonlinear states receive exactly the same held IMU input.
    X    = se23_propagate(X,   omega,a,g,dt);
    Xhat = se23_propagate(Xhat,omega,a,g,dt);
end

dxi = xi_exact_hist-xi_ll_hist;
dxM = xM_exact_hist-xM_lin_hist;

out.xi_exact = xi_exact_hist;
out.xi_ll    = xi_ll_hist;
out.xM_exact = xM_exact_hist;
out.xM_lin   = xM_lin_hist;

out.inv_att_mismatch  = vecnorm(dxi(1:3,:),2,1);
out.inv_vel_mismatch  = vecnorm(dxi(4:6,:),2,1);
out.inv_pos_mismatch  = vecnorm(dxi(7:9,:),2,1);
out.mekf_att_mismatch = vecnorm(dxM(1:3,:),2,1);
out.mekf_vel_mismatch = vecnorm(dxM(4:6,:),2,1);
out.mekf_pos_mismatch = vecnorm(dxM(7:9,:),2,1);
out.max_inv_component_error = max(abs(dxi),[],'all');
end

function xM = mekf_error_coordinates(X,Xhat)
%MEKF_ERROR_COORDINATES Exact coordinates for the tutorial convention:
% Rtilde=Rhat*R', vtilde=v-vhat, ptilde=p-phat.

[R,v,p] = se23_extract(X);
[Rhat,vhat,phat] = se23_extract(Xhat);
theta = so3_log_local(Rhat*R');
xM = [theta;v-vhat;p-phat];
end

function Xinv = se23_inv_local(X)
%SE23_INV_LOCAL Inverse of Gamma(R,v,p) on SE_2(3).
[R,v,p] = se23_extract(X);
Xinv = Gamma(R',-R'*v,-R'*p);
end

function xi = se23_log_local(X)
%SE23_LOG_LOCAL log(Gamma(R,v,p))^vee=[phi;J_l(phi)^(-1)v;J_l(phi)^(-1)p].
[R,v,p] = se23_extract(X);
phi = so3_log_local(R);
Jinv = so3_left_jacobian_inv_local(phi);
xi = [phi;Jinv*v;Jinv*p];
end

function phi = so3_log_local(R)
%SO3_LOG_LOCAL Numerically robust principal logarithm on SO(3).
c = min(1,max(-1,(trace(R)-1)/2));
theta = acos(c);

if theta < 1e-8
    phi = 0.5*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
    return;
end

% Principal-axis extraction near pi.
if abs(pi-theta) < 1e-5
    [V,D] = eig((R+eye(3))/2);
    [~,idx] = max(real(diag(D)));
    axis = real(V(:,idx));
    axis = axis/norm(axis);
    phi = theta*axis;
    return;
end

phi = theta/(2*sin(theta))* ...
      [R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)];
end

function Jinv = so3_left_jacobian_inv_local(phi)
%SO3_LEFT_JACOBIAN_INV_LOCAL Inverse SO(3) left Jacobian.
theta = norm(phi);
W = so3_hat(phi);

if theta < 1e-8
    Jinv = eye(3)-0.5*W+(1/12)*(W*W);
else
    c = 1/theta^2-(1+cos(theta))/(2*theta*sin(theta));
    Jinv = eye(3)-0.5*W+c*(W*W);
end
end

function format_panel(ax,T,fontSize,axesColor)
%FORMAT_PANEL Publication styling for the mismatch axes.
set(ax,'FontSize',fontSize,'LineWidth',1.1,'Color',axesColor, ...
       'GridAlpha',0.18,'MinorGridAlpha',0.12,'Layer','top');
xlim(ax,[0,T]);
end

function bound = one_significant_digit_upper_bound(value)
%ONE_SIGNIFICANT_DIGIT_UPPER_BOUND Conservative, readable annotation bound.
if ~(isfinite(value) && value > 0)
    bound = eps;
    return;
end
exponent = floor(log10(value));
% floor(scaled)+1 makes the displayed ``<'' bound strictly conservative,
% including the unlikely case that value is exactly a power-of-ten multiple.
mantissa = floor(value/10^exponent)+1;
if mantissa == 10
    mantissa = 1;
    exponent = exponent+1;
end
bound = mantissa*10^exponent;
end

function str = latex_precision_note(bound,unit)
%LATEX_PRECISION_NOTE Format a compact automatic numerical-precision note.
exponent = floor(log10(bound));
mantissa = round(bound/10^exponent);
str = sprintf('$\\mathrm{InvEKF}<%d\\times10^{%d}\\;\\mathrm{%s}$', ...
              mantissa,exponent,unit);
end

function print_case_summary(name,r)
%PRINT_CASE_SUMMARY Compact quantitative summary.
fprintf('\n%s motion:\n',name);
fprintf('  max |xi_exact-xi_LL| component : %.3e\n', ...
    r.max_inv_component_error);
fprintf('  final InvEKF velocity mismatch : %.3e m/s\n', ...
    r.inv_vel_mismatch(end));
fprintf('  final InvEKF position mismatch : %.3e m\n', ...
    r.inv_pos_mismatch(end));
fprintf('  final MEKF velocity mismatch   : %.3f m/s\n', ...
    r.mekf_vel_mismatch(end));
fprintf('  final MEKF position mismatch   : %.3f m\n', ...
    r.mekf_pos_mismatch(end));
end
