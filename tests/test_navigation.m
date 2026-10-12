function tests = test_navigation
%TEST_NAVIGATION Numerical identities for the tutorial's navigation tools.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
testCase.TestData.randomState = rng;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Tools'));
rng(2718,'twister');
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
rng(testCase.TestData.randomState);
end

function testHatOperator(testCase)
for k = 1:20
    w = randn(3,1); v = randn(3,1);
    W = so3_hat(w);
    verifyEqual(testCase,W*v,cross(w,v),'AbsTol',1e-14);
    verifyEqual(testCase,W+W',zeros(3),'AbsTol',1e-14);
end
end

function testExponentialAndIntegrals(testCase)
% A block matrix exponential independently evaluates R, J_l, and Q_l.
angles = unique([0, logspace(-12,1,160), ...
    1e-3*(1+[-1e-6,0,1e-6]), 0.1*(1+[-1e-6,0,1e-6]), pi, 2*pi]);
axes = [eye(3), [1;2;-3]/sqrt(14)];
maxError = zeros(1,3);
for axis = axes
    for theta = angles
        w = theta*axis;
        [Rref,Jref,Qref] = exponentialIntegrals(w);
        R = so3_exp(w);
        maxError = max(maxError,[norm(R-Rref,'fro'), ...
            norm(so3_left_jacobian(w)-Jref,'fro'), ...
            norm(so3_left_Q(w)-Qref,'fro')]);
        verifyEqual(testCase,R'*R,eye(3),'AbsTol',1e-13);
        verifyEqual(testCase,det(R),1,'AbsTol',1e-13);
    end
end
verifyLessThan(testCase,maxError,[1e-13,1e-12,1e-12]);
verifyEqual(testCase,so3_left_jacobian(zeros(3,1)),eye(3));
verifyEqual(testCase,so3_left_Q(zeros(3,1)),0.5*eye(3));
end

function testGroupExponential(testCase)
for k = 1:20
    xi = randn(9,1);
    Xi = zeros(5);
    Xi(1:3,1:3) = so3_hat(xi(1:3));
    Xi(1:3,4:5) = [xi(4:6),xi(7:9)];
    verifyEqual(testCase,se23_exp(xi),expm(Xi),'AbsTol',1e-13);
end
end

function testStateEmbedding(testCase)
R = so3_exp([0.2;-0.3;0.4]); v = [1;-2;3]; p = [4;2;-1];
X = Gamma(R,v,p);
[Rout,vout,pout] = se23_extract(X);
verifyEqual(testCase,Rout,R);
verifyEqual(testCase,vout,v);
verifyEqual(testCase,pout,p);
verifyEqual(testCase,X(4:5,:),[zeros(2,3),eye(2)]);
end

function testZeroRateAndZeroTime(testCase)
R = so3_exp([0.2;-0.3;0.4]); v = [1;-2;3]; p = [4;2;-1];
X = Gamma(R,v,p); a = [1.2;0.3;-9.6]; g = [0;0;9.81]; dt = 0.25;
expected = Gamma(R,v+(g+R*a)*dt,p+v*dt+0.5*(g+R*a)*dt^2);
verifyEqual(testCase,se23_propagate(X,zeros(3,1),a,g,dt),expected,'AbsTol',1e-14);
verifyEqual(testCase,se23_propagate_classical(X,zeros(3,1),a,g,dt),expected,'AbsTol',1e-14);
verifyEqual(testCase,se23_propagate(X,[1;2;3],a,g,0),X);
verifyEqual(testCase,se23_propagate_classical(X,[1;2;3],a,g,0),X);
end

function testConstantInputPropagation(testCase)
g = [0;0;9.81];
for k = 1:20
    R = so3_exp(randn(3,1)); v = randn(3,1); p = randn(3,1);
    X = Gamma(R,v,p); omega = randn(3,1); a = randn(3,1);
    dt = 0.03+0.3*rand;
    [E,J,Q] = exponentialIntegrals(omega*dt);
    expected = Gamma(R*E,v+g*dt+R*J*a*dt, ...
        p+v*dt+0.5*g*dt^2+R*Q*a*dt^2);
    Y = se23_propagate(X,omega,a,g,dt);
    verifyEqual(testCase,Y,expected,'AbsTol',1e-12);
    half = se23_propagate(X,omega,a,g,dt/2);
    half = se23_propagate(half,omega,a,g,dt/2);
    verifyEqual(testCase,Y,half,'AbsTol',1e-12);
end
end

function testInvariantErrorDynamics(testCase)
g = [0;0;9.81];
for k = 1:20
    X = Gamma(so3_exp(randn(3,1)),randn(3,1),randn(3,1));
    xi = randn(9,1); Xhat = se23_exp(xi)*X;
    omega = randn(3,1); a = randn(3,1); dt = 0.03+0.3*rand;
    Y = se23_propagate(X,omega,a,g,dt);
    Yhat = se23_propagate(Xhat,omega,a,g,dt);
    verifyEqual(testCase,Yhat/Y,se23_exp(Ak_RIEKF(g,dt)*xi),'AbsTol',1e-12);
end
end

function testTransitionMatrix(testCase)
g = [0;0;9.81];
A = zeros(9); A(4:6,1:3) = so3_hat(g); A(7:9,4:6) = eye(3);
for dt = [0,1e-5,0.01,0.1,1]
    verifyEqual(testCase,Ak_RIEKF(g,dt),expm(A*dt),'AbsTol',1e-13);
    verifyEqual(testCase,Ak_RIEKF(g,dt/2)^2,Ak_RIEKF(g,dt),'AbsTol',1e-13);
end
end

function testNoiseInjectionAdjoint(testCase)
% Differentiate group conjugation to check the continuous noise map.
X = Gamma(so3_exp([0.4;-0.3;0.2]),[1;-2;3],[4;2;-1]);
h = 1e-6; J = zeros(9,6);
for k = 1:6
    xi = zeros(9,1); xi(k) = h;
    D = (X*se23_exp(xi)/X-X*se23_exp(-xi)/X)/(2*h);
    J(:,k) = [D(3,2);D(1,3);D(2,1);D(1:3,4);D(1:3,5)];
end
dt = 0.1;
verifyEqual(testCase,Lk_RIEKF(X,[0;0;9.81],dt),dt*J,'AbsTol',1e-9);
end

function testCircularImu(testCase)
V = 5; r = 12.5; g = [0;0;9.81]; h = 1e-5;
for t = [0,1.3,7,20]
    [R,v,~] = circular_ground_truth(t,V,r);
    [Rp,vp,pp] = circular_ground_truth(t+h,V,r);
    [Rm,vm,pm] = circular_ground_truth(t-h,V,r);
    [omega,a] = imu_circular_motion(t,V,r,g);
    verifyEqual(testCase,(Rp-Rm)/(2*h),R*so3_hat(omega),'AbsTol',1e-9);
    verifyEqual(testCase,(vp-vm)/(2*h),R*a+g,'AbsTol',1e-9);
    verifyEqual(testCase,(pp-pm)/(2*h),v,'AbsTol',1e-9);
end
end

function testCircularTrajectory(testCase)
V = 5; r = 12.5; g = [0;0;9.81]; T = 20;
dtValues = [0.01,0.05,0.1]; classicalError = zeros(size(dtValues));
for j = 1:numel(dtValues)
    dt = dtValues(j);
    [R,v,p] = circular_ground_truth(0,V,r); X = Gamma(R,v,p); Xc = X;
    for k = 0:round(T/dt)-1
        [omega,a] = imu_circular_motion(k*dt,V,r,g);
        X = se23_propagate(X,omega,a,g,dt);
        Xc = se23_propagate_classical(Xc,omega,a,g,dt);
    end
    [R,v,p] = circular_ground_truth(T,V,r);
    verifyEqual(testCase,X,Gamma(R,v,p),'AbsTol',1e-9);
    classicalError(j) = norm(Xc(1:3,5)-p);
end
verifyGreaterThan(testCase,classicalError(1),0.01);
verifyGreaterThan(testCase,diff(classicalError),zeros(1,2));
end

function [R,J,Q] = exponentialIntegrals(w)
B = zeros(9);
B(1:3,1:3) = so3_hat(w);
B(1:3,4:6) = eye(3);
B(4:6,7:9) = eye(3);
E = expm(B);
R = E(1:3,1:3); J = E(1:3,4:6); Q = E(1:3,7:9);
end
