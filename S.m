function [ f] = S( u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[nx,ny]=size(u);

ax = 0; 
bx = 1;
ay = 0;
by = 1;
h = (bx-ax)/(nx-1);
h4 = h^4;
ii = 1:nx; x = ax + (ii-1)*h;
jj = 1:ny; y = ay + (jj-1)*h;
[X,Y] = ndgrid(x,y);
f = ffun( X, Y );

for j=2:ny-1
  for i=2:nx-1
    a1=(u(i+1,j)+u(i-1,j))/2;
    a2=(u(i,j+1)+u(i,j-1))/2;
    a3=(u(i+1,j+1)+u(i-1,j-1))/2;
    a4=(u(i-1,j+1)+u(i+1,j-1))/2;
    f(i,j) = (4*(a1-u(i,j))*(a2-u(i,j))-0.25*(a3-a4)^2)/h4;
  end
end



end

