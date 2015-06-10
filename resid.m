function r = resid( f, u, h )
%
% function r = resid( f, u, h )
%
% Computes the residual for the MOnge-Ampere problem
%
% Input
%   f      right-hand side
%   u      current approximate solution
%   h      mesh spacing

[nx,ny] = size( f );
r = zeros(nx,ny);
h4 = h^4;
for j=2:ny-1
  for i=2:nx-1
    a1=(u(i+1,j)+u(i-1,j))/2;
    a2=(u(i,j+1)+u(i,j-1))/2;
    a3=(u(i+1,j+1)+u(i-1,j-1))/2;
    a4=(u(i-1,j+1)+u(i+1,j-1))/2;
    r(i,j) = f(i,j) - (4*(a1-u(i,j))*(a2-u(i,j))-0.25*(a3-a4)^2)/h4;
  end
end
