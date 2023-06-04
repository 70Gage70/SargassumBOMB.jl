function [x, y] = sph2xy(lambda, lambda0, theta, theta0)
%SPH2XY Spherical to curvilinear spherical.
% [X Y] = SPH2SPH(LAMBDA,LAMBDA0,THETA,THETA0) where X,Y are in meters
% and LAMBDA's ,THETA's are in degrees.

R = 6371 * 1e3;
deg2rad = pi/180;
x = R * (lambda - lambda0)*deg2rad * cos(theta0*deg2rad);
y = R * (theta - theta0)*deg2rad;
