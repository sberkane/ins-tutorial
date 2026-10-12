function result = discretization(outputDir)
%DISCRETIZATION Exact and classical propagation for circular motion.

if nargin < 1, outputDir = []; end
[outputDir, cleanup] = experiment_setup(outputDir); %#ok<ASGLU>

% --------- GLOBAL FIGURE FONT SETTINGS (LaTeX-friendly) ---------
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

FONT_AXES   = 20;   % tick labels
FONT_LABEL  = 20;   % x/y labels
FONT_LEGEND = 20;
FONT_TITLE  = 22;
LINE_WIDTH  = 2;
bg_color = [1.0, 0.9, 0.9];   % approx red!10

% --- Scenario (circular motion) ---
T = 20;                     % seconds
g = [0;0;9.81];              % gravity in world (down)

V = 5.0;                     % m/s
w = 0.4;                     % rad/s
r = V/w;

% Initial conditions
R0 = eye(3);
v0 = [V;0;0];
p0 = [0;0;0];

% Sampling intervals
dt_list = [0.01 0.05 0.1];
exact = cell(size(dt_list)); classical = exact;

% Ground truth sampled finely (independent of dt)
dt_gt = 1e-3;
t_gt  = 0:dt_gt:T;
p_gt  = zeros(3,numel(t_gt));
for k = 1:numel(t_gt)
    [~,~,p_gt(:,k)] = circular_ground_truth(t_gt(k), V, r);
end

% --- Figure layout (two subfigures) ---
figure('Color', bg_color,'Units','pixels','Position',[100 100 1200 650]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Left: Exact SE2(3)

nexttile; hold on; grid on; axis equal;
ax = gca;
ax.FontSize = FONT_AXES;
ax.LineWidth = 1.2;

plot(p_gt(1,:), p_gt(2,:), 'k:', 'LineWidth', 3);  % Ground truth

for i = 1:numel(dt_list)
    dt = dt_list(i);
    p_ex = run_propagation(dt, T, R0, v0, p0, V, r, g, "exact");
    exact{i} = p_ex;
    plot(p_ex(1,:), p_ex(2,:), 'LineWidth', LINE_WIDTH);
end

xlabel('$x$ (m)','FontSize',FONT_LABEL);
ylabel('$y$ (m)','FontSize',FONT_LABEL);
lg = legend_strings(dt_list, "Exact");
legend(['Ground truth', lg], ...
       'Location','best', ...
       'FontSize',FONT_LEGEND);

% Right: Classical

nexttile; hold on; grid on; axis equal;
ax = gca;
ax.FontSize = FONT_AXES;
ax.LineWidth = 1.2;

plot(p_gt(1,:), p_gt(2,:), 'k:', 'LineWidth', 3);  % Ground truth

for i = 1:numel(dt_list)
    dt = dt_list(i);
    p_cl = run_propagation(dt, T, R0, v0, p0, V, r, g, "classical");
    classical{i} = p_cl;
    plot(p_cl(1,:), p_cl(2,:), 'LineWidth', LINE_WIDTH);
end

xlabel('$x$ (m)','FontSize',FONT_LABEL);
ylabel('$y$ (m)','FontSize',FONT_LABEL);

lg = legend_strings(dt_list, "Classical");
legend(['Ground truth', lg], ...
       'Location','best', ...
       'FontSize',FONT_LEGEND);

title(tl,'Exact Lie-group vs Classical discretization', ...
      'FontSize',FONT_TITLE,'FontWeight','bold');

axs = findall(gcf,'Type','axes');
linkaxes(axs,'xy');
% --- Axis padding (fraction of span) ---
pad_frac = 0.08;   % Horizontal padding

% Get current limits (now identical across axes)
xl = axs(1).XLim;
yl = axs(1).YLim;
dx = xl(2) - xl(1);
new_xl = xl + [-pad_frac, pad_frac]*dx;
set(axs,'XLim',new_xl,'YLim',yl);

exportgraphics(gcf,fullfile(outputDir,'exact_vs_classical.pdf'), ...
    'ContentType','vector','Padding',12);

result = struct('dt',dt_list,'T',T,'speed',V,'radius',r, ...
    'exact',{exact},'classical',{classical});
save(fullfile(outputDir,'discretization.mat'),'result');
end

function [p_out] = run_propagation(dt, T, R0, v0, p0, V, r, g, mode)
% Simulate propagation with a given dt, returning the estimated position trace.

t = 0:dt:T;
N = numel(t);

R = R0; v = v0; p = p0;
p_out = zeros(3,N);
p_out(:,1) = p;

for k = 1:N-1
    % Ideal IMU at time t(k)
    [omega, a] = imu_circular_motion(t(k), V, r, g);

    X = Gamma(R, v, p);

    switch mode
        case "exact"
            Xnew = se23_propagate(X, omega, a, g, dt);
        case "classical"
            Xnew = se23_propagate_classical(X, omega, a, g, dt);
        otherwise
            error('Unknown mode.');
    end

    [R, v, p] = se23_extract(Xnew);
    p_out(:,k+1) = p;
end
end

function lg = legend_strings(dt_list, prefix)
% Builds legend entries like "Exact dt=0.01s"
lg = strings(1,numel(dt_list));
for i = 1:numel(dt_list)
    lg(i) = sprintf('%s  $\\Delta t = %.3g\\,\\mathrm{s}$', prefix, dt_list(i));
end
lg = cellstr(lg);
end
