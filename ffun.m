function f = ffun( X, Y )
% FFUN  Forcing  f  of the Monge-Ampere problem
%
% Usage:  f = ffun( X, Y )
%
%     Input:
%         X  2D array of  x  coordinates of rectangular grid
%         Y  2D array of  y  coordinates of rectangular grid
%     Output:
%         f  contains values of the forcing on the grid

% evaluate the solution


f = (1.+X.^2+Y.^2).*exp(X.^2+Y.^2);
%f=2./(2-X.^2-Y.^2).^2;
%f=max(1-0.2./sqrt((X-0.5).^2+(Y-0.5).^2),0);
%f=1./sqrt(X.^2+Y.^2);

% h=X(2,1)-X(1,1);
% nx=1/h+1;
% ny=nx;
% epsilon=h;
% f=zeros(size(X));
% 
% for j=1:ny
% for i=1:nx
% if abs((X(i,1)-0.5)/epsilon)<1 && abs((X(i,1)-0.5)/epsilon)<1
% f(i,j)= -(abs((X(i,1)-0.5)/epsilon)+1)/epsilon;
% %f(i,j)= f(i,j)+(-(abs((X(i,1)-0.5)/epsilon)+1)/epsilon);
% else
% f(i,j)=0;
% end
% end
% end