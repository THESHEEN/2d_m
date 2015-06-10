function uc = injf2c( u )
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
nxc=(nx-1)/2+1;
nyc=(ny-1)/2+1;
uc=zeros(nxc,nyc);
uc(1,:)=u(1,1:2:ny);
uc(nxc,:)=u(nx,1:2:ny);
uc(:,1)=u(1:2:nx,1);
uc(:,nyc)=u(1:2:nx,ny);

for i=2:nxc-1
    for j=2:nyc-1
        uc(i,j)=1/4*u(2*i-1,2*j-1)+1/8*(u(2*i,2*j-1)+u(2*i-2,2*j-1)+u(2*i-1,2*j)+u(2*i-1,2*j-2))+1/16*(u(2*i,2*j)+u(2*i-2,2*j-2)+u(2*i,2*j-2)+u(2*i-2,2*j));
    end
end

%
% 2D full-weighting restriction operator
% 
% [nx,ny]=size()
% M = length(vm)-1;  
% Mym = M/2; ym = zeros(Mym+1);
% 
% ind = 3:2:M-1;  
% ym(2:Mym,2:Mym) = ( 4*vm(ind,ind) + ...
%                     2*(vm(ind-1,ind)+vm(ind+1,ind)+ ...
%                        vm(ind,ind-1)+vm(ind,ind+1)) + ...
%                       (vm(ind-1,ind-1)+vm(ind-1,ind+1)+ ...
%                        vm(ind+1,ind-1)+vm(ind+1,ind+1)) )/16;
% put ym in y 
% 
% y = ym(:);

