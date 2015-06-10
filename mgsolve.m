function [r,e] = mgsolve(q, N, etol, ncycles )
% MGSOLVE:  solves the Monge-Ampere problem by a multigrid method
%
% Usage:  [rnorm,enorm,sweeps] = mgsolve( N, rtol, ncycles );
%
% Input:
%    N        number of grid intervals in  x  and  y [default:  32]
%    rtol     relative tolerance for residual norm   [default:  1.0e-04]
%    ncycles  number of multigrid cycles to execute  [default:  10]
% Output:
%    rnorm    vector of residual norms after each cycle
%    enorm    vector of error    norms after each cycle
%    sweeps   vector of corresponding numbers of sweeps
%  (this code also plots  rnorm  and  enorm  vs.  sweeps)
%
% Uses:  UFUN, FFUN, RESID, MGV (which uses GSRELAX, INJF2C, CTOF)

% process the input arguments (set defaults if not supplied)
if (~exist('q','var')) || (isempty(q)), q = 0 ; end % ENO interpolaiton
if (~exist('N','var')) || (isempty(N)), N = 128 ; end
if (~exist('etol','var')) || (isempty(etol)), etol = 1.0e-8; end
if (~exist('ncycles','var')) || (isempty(ncycles)), ncycles = 100; end
clear rnorm enorm sweeps % just in case they were set before


n1 = 5; % number of sweeps on downward branch
n2 = 5; % number of sweeps on upward   branch
n3 = 5; % number of sweeps on coarsest grid


% set up the grid

ax = 0; bx = 1;
ay = 0; by = 1;
nx = N+1; h = (bx-ax)/(nx-1);
ny = nx; by = ay + (ny-1)*h;
ii = 1:nx; x = ax + (ii-1)*h;
jj = 1:ny; y = ay + (jj-1)*h;
[X,Y] = ndgrid(x,y);	% like meshgrid but x index is first, y is second

% set up the arrays for solution, forcing, and analytical solution

u = zeros(nx,ny);
f = zeros(nx,ny);
uexact = zeros(nx,ny);

% set up the analytical solution

uexact = ufun( X, Y );

% set up the forcing

f = ffun( X, Y );

% set up the boundary values

u(1:nx, 1) = uexact(1:nx, 1);
u(1:nx,ny) = uexact(1:nx,ny);

u( 1,1:ny) = uexact( 1,1:ny);
u(nx,1:ny) = uexact(nx,1:ny);
%u=uexact;
% solve by the multigrid method:  V-cycles


sweeps(1) = 0;
rnorm(1)  = norm( resid(f,u,h) )*h;
enorm(1)  = norm( u - uexact )*h;

for icycle=2:ncycles+1 % indexed from 2 since Matlab can't index arrays from 0
    
        u = mgv( q, f, u, h, n1, n2, n3 );
        sweeps(icycle) = sweeps(icycle-1) + (4*(n1+n2)/3);
        rnorm(icycle)  = norm( resid(f,u,h) )*h;
        %rnorm(icycle)  = norm(resid(f,u,h),inf );
        enorm(icycle)  = norm(u(:)-uexact(:),inf );
	cfacts = (rnorm(icycle)/rnorm(icycle-1))^(1/(n1+n2));
    
	disp(['cycle ', num2str(icycle-1),...
              ' error norm = ', num2str(enorm(icycle)),...
		' convergence/sweep = ', num2str(cfacts)]);
	%if abs(rnorm(icycle)-rnorm(icycle-1))<=etol, break, end; % quit loop if converged
    if norm(enorm(icycle)-enorm(icycle-1),inf)<=etol, break, end; % quit loop if converged
end
r=rnorm(icycle);
e=enorm(icycle);


% plot residual vs. number of sweeps
subplot(1,2,1);
semilogy(sweeps,rnorm,'r*-');
xlabel('Relaxation sweeps');
ylabel('Residual Norm');
title('Multigrid Method');

% plot solution error vs. number of sweeps
subplot(1,2,2);
semilogy(sweeps,enorm,'r*-');
xlabel('Relaxation sweeps');
ylabel('Error Norm');
title('Multigrid Method');

