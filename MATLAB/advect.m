function [xt, yt] = advect(xuv, yuv, tuv, u, v, x0, y0, tspan, tout)
%ADVECTION Plain advection.
% [XT, YT] = ADVECT(X, Y, T, U, V, X0, Y0, TSPAN, TOUT)

% F.J. Beron-Vera, 2021.
% Revised 2021/08/01 12:36:31.

if nargin < 9
   tout = tspan;
end

% create u,v interpolant
[xuv, yuv, tuv] = meshgrid(xuv, yuv, tuv);
P = [2 1 3]; % hay que permutar
xuv = permute(xuv, P);
yuv = permute(yuv, P);
tuv = permute(tuv, P);
u = permute(u, P);
v = permute(v, P);
Fu = griddedInterpolant(xuv, yuv, tuv, u);
Fv = griddedInterpolant(xuv, yuv, tuv, v);
clear xuv yuv tuv u v

% advection parameters
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% t0
xy0 = [x0(:); y0(:)];

% t
xy = ode4( ...
   @uv, ...
   tspan, xy0, ...
   Fu, Fv);
% [~, xy] = ode45( ...
%    @uv, ...
%    tspan, xy0, ...
%    options, ...
%    Fu, Fv);
xt = xy(:,1:end/2);
yt = xy(:,end/2+1:end);
xt = xt.';
yt = yt.';

if nargin > 9 % output
   Itout = findnearest(tspan, tout);
   xt = xt(:,Itout);
   yt = yt(:,Itout);
end

disp('MATLAB::Done!')

%--------------------------------------------------------------------------
function out = uv(t, xy, Fu, Fv)

n = length(xy)/2;
out = zeros(2*n,1);

x = xy(1:n);
y = xy(n+1:2*n);
t = repmat(t, [n 1]);

out(1:n) = Fu(x, y, t);
out(n+1:2*n) = Fv(x, y, t);