function [lambda, theta]=xy2sph(x, lambda0, y, theta0)
%XY2SPH Curvilinear spherical to spherical.
% [LAMBDA THETA] = XY2SHP(X,LAMBDA0,Y,THETA0) where X,Y are in meters
% and LAMBDA0,THETA0 are in degrees.

R = 6371 * 1e3;
deg2rad = pi/180;
lambda = lambda0 + x/(R*cos(theta0*deg2rad)) / deg2rad;
theta = theta0 + y/R / deg2rad;
