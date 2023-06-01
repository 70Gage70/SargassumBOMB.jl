function trajectories_ode4_mrs_short

%parametros
lat0 = 10;
lon0 = -90;
delta = 1.25; %3; %2.5;
a = 1*1e-4; %10cm de planta; 1*1e-4% flotador del sargaso 4mm en promedio (rango 2-7mm); 4 * 1e-6; %poner en  [km]
k = 20; %elastico [1/d]
n = 25;

[alfa,tau,R,f]=parametros(delta,a,lat0);

load viento_2021.mat
lonw=Lon;
latw=Lat;
tw=t;
uw=u; %km/d
vw=v;

[xa, ya] = sph2xy(lonw, lon0, latw, lat0);

xa=double(xa/1000);
ya=double(ya/1000);

clear Lon Lat t u v

load('merged-2021-IAS.mat')

% (u,v) [km/d]

[x, y] = sph2xy(lon, lon0, lat, lat0);
x=double(x/1000); %[km]
y=double(y/1000);

xdat=x;
ydat=y;
tdat=tw;

[xdat,ydat,tdat] = meshgrid(xdat,ydat,tdat);

clear lon lat

vx = interp3(x,y,t,u,xdat,ydat,tdat,'nearest');
vy = interp3(x,y,t,v,xdat,ydat,tdat,'nearest');

vax = interp3(xa,ya,tw,uw,xdat,ydat,tdat,'nearest');
vay = interp3(xa,ya,tw,vw,xdat,ydat,tdat,'nearest');

clear x y t u v xa ya tw uw vw

%zeros for current/wind on land
I=find(isnan(vx));
vax(I)=0;
vay(I)=0;
vx(I)=0;
vy(I)=0;


% Dv/Dt, Dva/Dt
[vx_x, vx_y, vx_t] = gradient(vx, xdat(1,:,1), ydat(:,1,1), tdat(1,1,:));
[vy_x, vy_y, vy_t] = gradient(vy, xdat(1,:,1), ydat(:,1,1), tdat(1,1,:));
DvxDt = vx_t + vx.*vx_x + vy.*vx_y;
DvyDt = vy_t + vx.*vy_x + vy.*vy_y;
[vax_x, vax_y, vax_t] = gradient(vax, xdat(1,:,1), ydat(:,1,1), tdat(1,1,:));
[vay_x, vay_y, vay_t] = gradient(vay, xdat(1,:,1), ydat(:,1,1), tdat(1,1,:));

% omega
omega = vy_x - vx_y;

% u & Du/Dt
ux = (1-alfa)*vx + alfa*vax;
uy = (1-alfa)*vy + alfa*vay;
ux_t = (1-alfa)*vx_t + alfa*vax_t;
uy_t = (1-alfa)*vy_t + alfa*vay_t;
ux_x = (1-alfa)*vx_x + alfa*vax_x;
ux_y = (1-alfa)*vx_y + alfa*vax_y;
uy_x = (1-alfa)*vy_x + alfa*vay_x;
uy_y = (1-alfa)*vy_y + alfa*vay_y;
DuxDt = ux_t + ux.*ux_x + uy.*ux_y;
DuyDt = uy_t + ux.*uy_x + uy.*uy_y;

xdat=squeeze(xdat(1,:,1))';
ydat=squeeze(ydat(:,1,1));
tdat=squeeze(tdat(1,1,:));


%x0,y0
load IC

%{
 [X0 Y0]=meshgrid(xdat,ydat);
%  [X0, Y0] = sph2xy(X0, lon0, Y0, lat0);

I=find(vx(:,:,1)==0 & vy(:,:,1)==0);
X0(I)=[];
Y0(I)=[];
X0 = X0(:); %*1e-3; % [km]
Y0 = Y0(:); %*1e-3; % [km]
I=find(X0<615 & Y0<597);
X0(I)=[];
Y0(I)=[];

I=find(X0<-232 & Y0<820);
X0(I)=[];
Y0(I)=[];

I=find(X0<1355 & Y0<-264);
X0(I)=[];
Y0(I)=[];

I=find(X0<1437 & Y0<-236);
X0(I)=[];
Y0(I)=[];
%load linea
I = inpolygon(X0, Y0, xb, yb); I = find(I>0);
X0(I)=[];
Y0(I)=[];
%}

% [X0, Y0] = sph2xy(-78, lon0,15, lat0);
% X0 = X0(:)*1e-3; % [km]
% Y0 = Y0(:)*1e-3; % [km]

t0 = datenum(2021,2,1);
T= 2;
length(X0)
load BOM_T150days_
load BOM_T150days

%N1=[];

for t0_= t0+15:15:t0+150

tz=t0_:.1:t0_+T;

X0_=[];
Y0_=[];
XT_=[];
YT_=[];

    for N= [N1 N2] %10:20:length(X0)

        disp('N =')
        N
        [x0, y0, edges] = net(X0(N), Y0(N), .5, .5, 2, n, 0);
        ell = sqrt(diff(x0(edges')).^2 + diff(y0(edges')).^2)';

        z0 = [x0(:);  y0(:)];

        z = ode4(@ebom_slow, ...
            tz, z0, ...
            vx, vy, DvxDt, DvyDt, omega, ux, uy, ...
            DuxDt, DuyDt, ...
            xdat, ydat, tdat, ...
            R, tau, f,k,ell,edges);

        x = z(:,1:end/2)';
        y = z(:,end/2+1:end)';

        x=nanmean(x);
        y=nanmean(y);

        [x, y] = xy2sph(x*1e3, lon0, y*1e3, lat0); % [deg]

        X0_=[X0_; x(:,1)];
        XT_=[XT_; x(:,end)];
        Y0_=[Y0_; y(:,1)];
        YT_=[YT_; y(:,end)];

    end

    eval(['save BOM_T2days_' num2str(t0_)  ' t0 X0_ Y0_ XT_ YT_ alfa tau a delta k n'])

end

disp('MATLAB::DONE!')

%-------------------------------------------------------------------------------------------------------------------
function [alfa,tau,R,f]=parametros(delta,a,lat0)

% parameters
f = 2 * 2*pi * sin(lat0*pi/180);        % [1/d]
rho = 1027 * 1e9;                             % [kg/km^3]
rhoa = 1.2 * 1e9;                             % [kg/km^3]
nu = 1e-6 * .24*.36;                          % [km^2/d]
nua = 1.5e-5 * .24*.36;                       % [km^2/d]
mu = nu*rho;
mua = nua*rhoa;
gamma = mua/mu;

 i=sqrt(-1);
    phi = (i*sqrt(1-(2./delta-1).^2)+2./delta-1).^(1/3);
    Phi = i*sqrt(3)/2*(1./phi-phi)-1./(2*phi)-phi/2+1;
    Phi = real(Phi);
    R = (1 - Phi/2)/(1 - Phi/6);
    Psi = 1/pi*acos(1-Phi)-1/pi*(1-Phi).*sqrt(1-(1-Phi).^2);

    alfa = gamma*Psi./(gamma*Psi+(1-Psi));
    tau = ((1 - Phi/6)) / ((1 + (gamma-1)*Psi)*delta^4) * a^2/(3*mu/rho);

%-------------------------------------------------------------------------------------------------------------------
function [x, y, edges] = net(x0, y0, dx, dy, dimension, gridsz, A)
% NET Net genereation.

% Philippe Miron 2019/11/01

% node and edge structure
% 4x4 nbpts representing sargassum's mat
%    x1 - x2 - x3
%    |    |    |
%    x5 - x6 - x7
%    |    |    |
%    x9 - x10 -x11
%    ^
%    |
% (x0,y0)

if nargin < 7
        A = 0;
end

% number of points
if dimension == 1 % chain
     nbpts = gridsz;
     x = ones(nbpts, 1)*x0 + (1:nbpts)'*dx + A*rand(nbpts,1);
     y = ones(nbpts, 1)*y0 + A*rand(nbpts,1);
elseif dimension == 2
     nbpts = gridsz*gridsz;
     x = ones(nbpts, 1)*x0 + repmat(0:gridsz-1,[1,gridsz])'*dx +  A*rand(nbpts,1);
     y = ones(nbpts, 1)*y0 + repelem(0:gridsz-1,gridsz)'*dy + A*rand(nbpts,1);
end

% edge and node structure
% the edges matrix contain all the links
% horizontal and vertical between nodes
% in the grid
% edges(1,:) = [n1, n2];
% edges(2,:) = [n2, n3];
% edges structure creation
edges = [];
if dimension == 1
     for i = 1:gridsz-1
         edges = [edges; [i, i+1]];
     end
elseif dimension == 2
     % horizontal connection
     for i = 1:gridsz
         for j = 1:gridsz-1
             ind = (i-1)*gridsz + j;
             edges = [edges; [ind,ind+1]];
         end
     end

     % vertical connection
     for i = 1:gridsz-1
         for j = 1:gridsz
             ind = (i-1)*gridsz + j;
             edges = [edges; [ind,ind+gridsz]];
         end
     end
end
%-------------------------------------------------------------------------------------------------------------------

function out = ebom_slow(t, in, vx, vy, DvxDt, DvyDt, omega, ux, uy, DuxDt, DuyDt, xdat, ydat, tdat, R, tau, f, k, ell, edges)

n = length(in)/2;
out = zeros(2*n,1);

x = in(1:n);
y = in(n+1:2*n);

nx = interp1(xdat, 1:length(xdat), x);
ny = interp1(ydat, 1:length(ydat), y);
nt = interp1(tdat, 1:length(tdat), t);
nt = repmat(nt, [n 1]);

vx = ba_interp3(vx, nx, ny, nt, 'cubic');
vy = ba_interp3(vy, nx, ny, nt, 'cubic');
DvxDt = ba_interp3(DvxDt, nx, ny, nt, 'cubic');
DvyDt = ba_interp3(DvyDt, nx, ny, nt, 'cubic');
omega = ba_interp3(omega, nx, ny, nt, 'cubic');
ux = ba_interp3(ux, nx, ny, nt, 'cubic');
uy = ba_interp3(uy, nx, ny, nt, 'cubic');
DuxDt = ba_interp3(DuxDt, nx, ny, nt, 'cubic');
DuyDt = ba_interp3(DuyDt, nx, ny, nt, 'cubic');

K = zeros(n);
for iedge = 1:size(edges,1)
    n1 = edges(iedge, 1);
    n2 = edges(iedge, 2);
    K(n1, n1) = K(n1, n1) - k;
    K(n1, n2) = K(n1, n2) + k;
    K(n2, n2) = K(n2, n2) - k;
    K(n2, n1) = K(n2, n1) + k;
end
L = zeros(n);
xij = [diff(x(edges'))' diff(y(edges'))'];
dij = vecnorm(xij')';
for iedge = 1:size(edges,1)
    n1 = edges(iedge,1);
    n2 = edges(iedge,2);
    Li = k * ell(iedge) / dij(iedge);
    L(n1, n1) = L(n1, n1) - Li;
    L(n1, n2) = L(n1, n2) + Li;
    L(n2, n2) = L(n2, n2) - Li;
    L(n2, n1) = L(n2, n1) + Li;
end
Fx = (K - L) * x;
Fy = (K - L) * y;

out(1:n)     = ux + tau * (R*DvxDt - R*(f + omega/3) .* vy - DuxDt + (f + R*omega/3).* uy + Fx);
out(n+1:2*n) = uy + tau * (R*DvyDt + R*(f + omega/3) .* vx - DuyDt - (f + R*omega/3).* ux + Fy);

