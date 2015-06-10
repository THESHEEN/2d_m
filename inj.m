function uc = inj( u )
%
% function uc = injf2c( u )
%
% Transfers a fine grid to a coarse grid by injection
%
% Input
%   u    fine-grid function
%
% Returns
%   uc   coarse-grid function
%

[nx,ny] = size( u );
uc = u(1:2:nx,1:2:ny);

