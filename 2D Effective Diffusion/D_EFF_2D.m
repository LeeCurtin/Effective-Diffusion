function [u1,p,t] = D_EFF_2D(tspan,gd,sf,ns)
%D_EFF_2D - Using a finite element method to evaluate the diffusion of
%matter out of/into squares of slow diffusion representing PLGA on a small scale.

tic

%CREATE MESH AND BOUNDARY
g = decsg(gd,sf,ns); %Creates g using the matrices produced by PDETOOLBOX
[p,e,t] = initmesh(g, 'hmax', 0.05);

%SET INITIAL CONDITIONS

[C,A,F]=assema(p,t,0,1,'0!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1!1'); 
u0=A\F;

% CREATE A PDE ENTITY FOR ONE PDE WITH ONE DEPENDENT VARIABLE
numberOfPDE = 1;
pb = pde(numberOfPDE);
pg = pdeGeometryFromEdges(g); % Creates a geometry entity

% BOUNDARY CONDITIONS
% Square = pdeBoundaryConditions(pg.Edges([1,2,3,4]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle in Middle
% Square = pdeBoundaryConditions(pg.Edges([1:5]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle on Bottom
% Circle = pdeBoundaryConditions(pg.Edges([5:8]),'q', [0,0,0;0,0,0;0,0,0], 'g', [0;0;0]); %Zero flux for Circle in circle
% pb.BoundaryConditions = [Square];

%SET COEFFICIENTS OF PDEs

%DIFFUSION COEFFICIENTS
D_water = 3600*50.e-5;
D_paste = 3600*67.e-7;

c_temp = sprintf('%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f!%1.3f',D_water,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste,D_paste);

c = char(c_temp);

%COEFFICIENTS OF A
a = ['0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00!0.00'];

%COEFFICIENTS OF F
% f = char('-0.5*u(1).*u(3)!-0.5*u(1).*u(3)','0.5*u(1).*u(3)!0.5*u(1).*u(3)','0.0!0.0');
f = char('0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0!0');

%COEFFICIENTS OF D
d = [1];

%SOLVE SYSTEM OF PDEs
u1 = parabolic(u0,tspan,pb,p,e,t,c,a,f,d);
toc

return




