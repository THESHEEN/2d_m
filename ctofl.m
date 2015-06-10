function [y] = ctofl( v )
%
% function u = ctof( uc )
%
% Transfers a coarse grid to a fine grid
%
% Input
%   uc    coarse-grid function
%
% Returns
%   u     fine-grid function

[nx,ny]=size(v);
M = nx-1;  
vm = reshape(v,M+1,M+1);
Mym = 2*M; ym = zeros(Mym+1);
nyf=Mym+1;nxf=nyf;

% % transfer by copying where the grids line up
% for j=1:ny
%   for i=1:nx
%     u(2*(i-1)+1,2*(j-1)+1) = uc(i,j);
%   end
% end

% ax = 0; 
% bx = 1;
% ay = 0;
% by = 1;
% h = (bx-ax)/Mym;
% ii = 1:nxf; x = ax + (ii-1)*h;
% jj = 1:nyf; y = ay + (jj-1)*h;
% [X,Y] = ndgrid(x,y);
% uexact = ufun( X, Y );

% ym(1:nxf, 1) = uexact(1:nxf, 1);
% ym(1:nxf,nyf) = uexact(1:nxf,nyf);
% ym( 1,1:nyf) = uexact( 1,1:nyf);
% ym(nxf,1:nyf) = uexact(nxf,1:nyf);

% ym center points 

ym(3:2:Mym-1,3:2:Mym-1) = vm(2:M,2:M);

% ym SE, NE, SW & NW points

ym(2:2:Mym,2:2:Mym) = ( vm(1:M,  1:M) + vm(2:M+1,  1:M) + ...
                        vm(1:M,2:M+1) + vm(2:M+1,2:M+1) )/4;
% ym N & S points                    

ym(3:2:Mym-1,2:2:Mym) = ( vm(2:M,1:M) + vm(2:M,2:M+1) )/2;

% ym E & W points
                    
ym(2:2:Mym,3:2:Mym-1) = ( vm(1:M,2:M) + vm(2:M+1,2:M) )/2;

% put ym in y 

y = ym;



% % linear interpolation in x
% for j=2:nyf-1
%   for i=2:2:nxf-1
%     u(i,j) = 1/4*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1));
%   end
% end
% 
% % linear interpolation in y
% for j=2:2:nyf-1
%   for i=2:nxf-1
%     u(i,j) = 1/4*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1));
%   end
% end
% 
