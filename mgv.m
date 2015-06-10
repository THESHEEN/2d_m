function u = mgv(q, f, u, h, n1, n2, n3, io )
% MGV  Does one V-cycle for the Monge-Ampere problem on a grid with mesh size  h
%
% Usage:  u = mgv( f, u, h, n1, n2, n3, io )
%
% Input
%    f    right-hand side
%    u    current approximate solution
%    h    mesh spacing
%    n1   number of sweeps (downward branch)
%    n2   number of sweeps (upward   branch)
%    n3   number of sweeps (coarsest grid)
%    io   [optional] use io>0 to print standard multigrid trace
%
% Returns
%    u    updated
%
% Uses:  GSRELAX, RESID, INJF2C, CTOF

if (~exist('io','var')) || (isempty(io)), io = 0; end

% set up the finest grid
if  size(f)~=size(u)
   disp('sizes of f and u must match');
   return
end
l = 1;
[nx(l),ny(l)] = size( f );
hval(l) = h;
eval(['f_' num2str(l) ' = f;']);
eval(['u_' num2str(l) ' = u;']);

% set up the coarse grids

while (mod(nx(l)-1,2)==0 && mod(ny(l)-1,2)==0 && nx(l)-1>2 && ny(l)-1>2 && l==1 )
   l = l+1;
   nx(l) = (nx(l-1)-1)/2 + 1;
   ny(l) = (ny(l-1)-1)/2 + 1;
   hval(l) = 2*hval(l-1);
   eval(['f_' num2str(l) ' = zeros(nx(l),ny(l));']);
   eval(['u_' num2str(l) ' = zeros(nx(l),ny(l));']);
end
nl = l;

% downward branch of the V-cycle
for l=1:nl-1
   ul = eval(['u_' num2str(l)]);
   fl = eval(['f_' num2str(l)]);
   for isweep=1:n1
      ul = gsrelax( fl, ul, hval(l) );
      if io>0
         rn = norm(resid( fl, ul, hval(l) ))*hval(l);
         disp(['level ',num2str(l),' residual norm = ', num2str(rn)]);
      end
   end

   r = resid( fl, ul, hval(l) );
   eval(['u_' num2str(l) ' = ul;']);
   eval(['u_' num2str(l+1) ' = injf2c( ul );']);

   eval(['f_' num2str(l+1) ' = injf2c( r )+ S(inj(ul));']);
  
end

% solve on the coarsest grid by many relaxation sweeps
l = nl;
ul = eval(['u_' num2str(l)]);
fl = eval(['f_' num2str(l)]);
for isweep=1:n3
   ul = gsrelax( fl, ul, hval(l) );
   if io>0
      rn = norm(resid( fl, ul, hval(l) ))*hval(l);
      disp(['level ',num2str(l),' residual norm =', num2str(rn)]);
   end
end
eval(['u_' num2str(l) ' = ul;']);

% upward branch of the V-cycle
for l=nl-1:-1:1
   uc = eval(['u_' num2str(l+1)]);
   ul = eval(['u_' num2str(l)]);
   if q==1
       ul = ul + ctofo( uc - inj(ul) );
   else
       ul = ul + ctofl( uc - inj(ul) );
   end    
   fl = eval(['f_' num2str(l)]);
   for isweep=1:n2
      ul = gsrelax( fl, ul, hval(l) );
      if io>0
         rn = norm(resid( fl, ul, hval(l) ))*hval(l);
         disp(['level ',num2str(l),' residual norm =', num2str(rn)]);
      end
   end
   eval(['u_' num2str(l) ' = ul;']);
end

% copy solution back to original grid
eval(['u = u_' num2str(1) ';']);
