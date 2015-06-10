function u = ufun( X, Y )
% UFUN  Analytical solution  u  of the Monge-Ampere problem
%
% Usage:  u = ufun( X, Y )
%
%     Input:
%         X  2D array of  x  coordinates of rectangular grid
%         Y  2D array of  y  coordinates of rectangular grid
%     Output:
%         u  contains values of the analytical solution on the grid
%
% See also:  MGSOLVE, FFUN


% evaluate the solution


u = exp((X.^2+Y.^2)./2);
%u=-sqrt(2-X.^2-Y.^2);
%u=0.5*max(sqrt((X-0.5).^2+(Y-0.5).^2)-0.2,0).^2;

%u = (2*sqrt(2)/3).*(X.^2+Y.^2).^0.75;

%u= 1.5.*min(X,0.5)-0.5.*max(X-0.5,0)%+1.5.*min(Y,0.5)-0.5*max(Y-0.5,0);