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

r=resid(f,u,h);
plot(abs(fft(r(20,:),N)),'b','Linewidth',2);
hold on

u=mgv(0,f,u,h,100,0,1);
r0=resid(f,u,h);
plot(abs(fft(r0(20,:),N)),'r','Linewidth',2);

u=mgv(1,f,u,h,100,0,1);
r1=resid(f,u,h);
plot(abs(fft(r1(20,:),N)),'c','Linewidth',2);