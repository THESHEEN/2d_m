N=32;
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

e=u-uexact;
plot(abs(fft(e(20,:),N)),'b','Linewidth',2);
hold on

% GS iteration
e=MA(1,N,191);
plot(abs(fft(e(20,:),N)),'r','Linewidth',2);

% Jacobi iteration
e=MA(0,N,191);
plot(abs(fft(e(20,:),N)),'g','Linewidth',2);

% ADI
e=MA(2,N,3);
plot(abs(fft(e(20,:),N)),'c','Linewidth',2);