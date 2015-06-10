function u = ctofo( uc )
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

[nx,ny] = size( uc );
nxf = 2*(nx-1)+1;
nyf = 2*(nx-1)+1;
u = zeros(nxf,nyf);

% transfer by copying where the grids line up
for j=1:ny
  for i=1:nx
    u(2*(i-1)+1,2*(j-1)+1) = uc(i,j);
  end
end

%ENO interpolation in x
for j=2:nyf-1
  for i=2:2:nxf-1
      if i==2
         u(i,j)=(-1/8)*u(i+3,j)+3/4*u(i+1,j)+3/8*u(i-1,j);
      elseif i==nxf-1
         u(i,j)=(-1/8)*u(i-3,j)+3/4*u(i-1,j)+3/8*u(i+1,j);
      elseif  u(i-3,j)-2*u(i-1,j)+u(i+1,j) <= u(i-1,j)-2*u(i+1,j)+u(i+3,j)
             u(i,j)=(-1/8)*u(i-3,j)+3/4*u(i-1,j)+3/8*u(i+1,j);
      else
            u(i,j)=(-1/8)*u(i+3,j)+3/4*u(i+1,j)+3/8*u(i-1,j);
      end
  end
end

%ENO interpolation in y
for j=2:2:nyf-1
  for i=2:nxf-1
      if j==2
          u(i,j)=(-1/8)*u(i,j+3)+3/4*u(i,j+1)+3/8*u(i,j-1);
      elseif j==nyf-1
         u(i,j)=(-1/8)*u(i,j-3)+3/4*u(i,j-1)+3/8*u(i,j+1);
      elseif abs(u(i,j-3)-2*u(i,j-1)+u(i,j+1))<= abs(u(i,j-1)-2*u(i,j+1)+u(i,j+3))
          u(i,j)=(-1/8)*u(i,j-3)+3/4*u(i,j-1)+3/8*u(i,j+1);
      else
          u(i,j)=(-1/8)*u(i,j+3)+3/4*u(i,j+1)+3/8*u(i,j-1);
      end
  end
end