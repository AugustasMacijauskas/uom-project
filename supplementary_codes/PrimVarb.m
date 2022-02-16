
%% Below is a Matlab code for the primitive variable formulation with a staggered grid.

%%  First try some tests to check the stability of time-stepping, the second-order accuracy, and the time to reach the steady state.
    
%%  Gather the results for different spatial resolutions $nx$ for the steady horizontal velocity at the mid-section $x=0.5$, and plot these results on top of one another. Include a result from the streamfunction vorticity formulation.

%%  The code calculates the force to second-order accuracy on the top plate, on the sides and on the bottom. 

%%   Experiment with different geometries of the cavity, say $ny=2*nx$. 
%      The optimal parameter for the SOR relaxation will change with the geometry, but this should not affect the calculation of the steady state.
%      The force on the bottom plate should decrease rapidly with deeper cavities, while the force on the top plate should decrease only slightly. On the other hand for shallower cavities, the force on the top plate increases as the height decreases.

%%   Experiment with different Reynolds numbers, $1/mu$, say $mu=0.2$, 0.5, and 0.05. The code reduces the time-step as $mu$ is increased. That will increase the run-time, although the steady state should be reached earlier. At smaller values of $mu$ the time-step may be limited by the CFL condition if the spatial resolution is not increased. The forces change very little in $mu>0.1$, i.e.~are effectively in the low Reynolds number limit. For higher Reynolds numbers, $mu<0.02$, the forces increase as the top boundary layer begins to thin. Note at these higher Reynolds number the steady state is achieved later, say by $t=10$. 



 %% primative variable on staggered grid, algorithm 3
 clear all %clear all variables
 clc % clear output screen

 %% Constants
   
 mu = 0.1;  %i.e. Re=10
 nx = 11+1; % +1 because index 1..nx not 0..nx
 ny = nx;  %different nx & ny gives non-square
 tfinal = 3.0;
 dtpic = 0.5;

 w = zeros(nx,ny); % form a matrix of zeros (1..nx, 1..ny) 
 p = w; % w = vorticty, p = streamfunction
 u = w; % u velocty, sim v
 v = w;
 ut = w; % ut = du/dt, sim vt
 vt = w;
 bcb = zeros(nx); % for BC on Bottom. Top, Left and Right
 bct = bcb;
 bcl = zeros(ny);
 bcr = bcl;
   
 %% Initialization
 h = 1/(nx-2); %h = dx (=dy); so x-size is 1, y-size is ny/nx
 t = 0;
 r = 2/(1+ 0.8/(nx-1)); %SOR parameter
 dt = 0.2*h*h/mu; 
 tpic = t + dtpic - 0.5*dt;
   
 for i= 1:nx
   for j=1:ny
     p(i,j) = 0;
     u(i,j) = 0;
     v(i,j) = 0;
   end
 end
   
 %% Main program starts here
 while(t<tfinal) %time-step loop starts here

 %% first part with no pressure, u then v   
 for i =2:nx-2      
   u(i,ny)= 2*sin(pi*(i-1)*h)*sin(pi*(i-1)*h) - u(i,ny-1);

   % no-slip BC using ghost point just outside
   u(i,1) = - u(i,2);
   u(i,ny)= 2*sin(pi*(i-1)*h)*sin(pi*(i-1)*h) - u(i,ny-1);
   % find du/dt in interior
   for j = 2:ny-1
     ut(i,j) = - u(i,j)*(u(i+1,j)-u(i-1,j))/(2*h) - (v(i,j)+v(i+1,j)+v(i,j-1)+v(i+1,j-1))*(u(i,j+1)-u(i,j-1))/(8*h) + mu*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/(h*h);
   end
 end
   
 for j = 2:ny-2
   % no-slip BC using ghost point just outside
   v(1,j) = - v(2,j);
   v(nx,j)= - v(nx-1,j);
   % find dv/dt in interior
   for i = 2:nx-1
	     vt(i,j) = - (u(i,j)+u(i,j+1)+u(i-1,j)+u(i-1,j+1))*(v(i+1,j)-v(i-1,j))/(8*h) - v(i,j)*(v(i,j+1)-v(i,j-1))/(2*h) + mu*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j))/(h*h);
   end
 end
   
 %% now time-step u and v in interior (no pressure)
 for i=2:nx-2
   for j=2:ny-1
     u(i,j)= u(i,j) + dt*ut(i,j);
   end
 end
 for i= 2:nx-1
   for j= 2:ny-2 
     v(i,j)= v(i,j) + dt*vt(i,j);
   end
 end

 %% BC for dp/dn = mu*d2u/dn2
 for i = 2:nx-1 
   bcb(i) = - mu*(-v(i,4)+4*v(i,3)-5*v(i,2)+2*v(i,1))/h;
   bct(i) = mu*(-v(i,ny-4)+4*v(i,ny-3)-5*v(i,ny-2)+2*v(i,ny-1))/h;
 end
 for j = 2:ny-1
   bcl(j) = - mu*(-u(4,j)+4*u(3,j)-5*u(2,j)+2*u(1,j))/h;
   bcr(j) = mu*(-u(nx-4,j)+4*u(nx-3,j)-5*u(nx-2,j)+2*u(nx-1,j))/h;
 end
 
 %% now time-step boundary u & v (no pressure) by above bc
 for i = 2:nx-1 
   v(i,1) = v(i,1) - dt*bcb(i)/h;
   v(i,ny-1) = v(i,ny-1) + dt*bct(i)/h;
 end
 for j = 2:ny-1
   u(1,j) = u(1,j) - dt*bcl(j)/h;
   u(nx-1,j) = u(nx-1,j) + dt*bcr(j)/h;
 end

 %% find (div u)/dt at interor points
 for i = 2:nx-1 
   for j = 2:ny-1 
     w(i,j) = (u(i,j) - u(i-1,j) + v(i,j) - v(i,j-1))/(dt*h);  
   end
 end 
 
 %% Poisson solver starts here -----------------
 %%    to solve nabla^2 p = - w  (= (div u)/dt)
 for it = 1:2*nx % iterations start here

   %% BC dp/dn=mu*d2un/dn2 - values found above in bcX
   for i = 2:nx-1
     p(i,1) = p(i,2) + bcb(i);
     p(i,ny)= p(i,ny-1) + bct(i);
   end              
   for j = 1:ny-1
     p(1,j)= p(2,j) + bcl(j);
     p(nx,j)= p(nx-1,j) + bcr(j);
   end              
           
   %% one iteration of p at interior points 
   for i = 2:nx-1
     for j = 2:ny-1
       p(i,j) =(1-r)*p(i,j) + ...
           r*0.25*(p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - w(i,j)*h*h );
     end
   end
 end    % iterations end here
           
 %% make average pressure zero
 q = 0;
 for i = 1:nx
   for j = 1:ny
     q = q + p(i,j);
    end  
  end
  q = -q/((nx)*(ny));
  for i = 1:nx
    for j = 1:ny
      p(i,j) = p(i,j) + q;
    end
  end
  %% Poisson solver Ends here -------------------

  %% second part of time-step - the pressure projection 
  %%    at interior points 
  for i =1:nx-1 % include left and right boundaries
    for j = 2:ny-1
      u(i,j) = u(i,j) - dt*(p(i+1,j) - p(i,j))/h;
    end
  end
  for j = 1:ny-1 % include top and bottom boundaries
    for i = 2:nx-1
      v(i,j) = v(i,j) - dt*(p(i,j+1) - p(i,j))/h;
    end
  end
       
  %% every dtpic print out force
  if t>tpic
    % Force calculation       
    % top, uses u = 0 at ends to save correction to int  
    q = 0;
    for i = 1:nx-1 
      q = q + (u(i,ny)-u(i,ny-1));          
    end     
    % bottom
    qb = 0;
    for i = 1:nx-1
      qb = qb - (u(i,2)-u(i,1));
    end      
    % sides       
    qs = 0;
    for j = 2:ny-1
      qs = qs + h*0.5*(p(1,j)+p(1,j) - p(nx-1,j)-p(nx,j))/mu;
    end
    formatSpecA='t= %8.4f; top= %8.4f;sides=%8.4f;bottom=%8.4f; -qs-qb=%8.4f\n ';
    fprintf(formatSpecA, t, q, qs, qb, -qs-qb);
    % End of Force calculation
    tpic = t + dtpic - 0.5*dt;
  end % endif t>tpic

  t = t + dt;
  end %% end of time-step, endwhile of t<tfinal
 
 %% velocity u on mid-section x=0.5
 figure(1)
 py = zeros(ny);
 pu = py;
 py(1) = 0;
 pu(1) = 0;
 py(ny) = h*(ny-2);
 pu(ny) = 1;
 for j = 2:ny-1
   py(j) = h*(j-1-0.5);
   pu(j) = u(round((nx-1)/2),j);
 end
 plot(pu,py,'-*');
 title(' u at mid-section x=0.5');
 xlabel('u');
 ylabel('y');

   
 %% Plot 'u' 
 figure(2)   
 tx = linspace (0, (nx-1)*h, nx)';
 ty = linspace (0, (ny-1-0.5)*h, ny)';
 [xx, yy] = meshgrid (tx, ty);
 mesh (tx, ty, u');
 xlabel ('x');
 ylabel ('y');
 zlabel ('u(x,y)');
 title ('u plot');
    
 %% Plot 'v'
 figure(3)
 tx = linspace (0, (nx-1-0.5)*h, nx)';
 ty = linspace (0, (ny-1)*h, ny)';
 [xx, yy] = meshgrid (tx, ty);
 mesh (tx, ty, v');
 xlabel ('x');
 ylabel ('y');
 zlabel ('v(x,y)');
 title ('v plot');
    
 %% Plot 'p'
 %% smooth corner points, not set elsewhere
 p(1,1) = p(2,1) + p(1,2) - p(2,2);
 p(nx,1) = p(nx-1,1) + p(nx,2) - p(nx-1,2);
 p(1,ny) = p(2,ny) + p(1,ny-1) - p(2,ny-1);
 p(nx,ny) = p(nx-1,ny) + p(nx,ny-1) - p(nx-1,ny-1);
 figure(4)
 tx = linspace (0, (nx-1-0.5)*h, nx)';
 ty = linspace (0, (ny-1-0.5)*h, ny)';
 [xx, yy] = meshgrid (tx, ty);
 mesh (tx, ty, p');
 xlabel ('x');
 ylabel ('y');
 zlabel ('p(x,y)');
 title ('p plot');

 %% Plot 'div u'
 %% find (div u)/dt at interor points
 for i = 2:nx-1 
   for j = 2:ny-1 
     w(i,j) = (u(i,j) - u(i-1,j) + v(i,j) - v(i,j-1))/h;  
   end
 end   
 figure(5)
 tx = linspace (0, (nx-1-0.5)*h, nx)';
 ty = linspace (0, (ny-1-0.5)*h, ny)';
 [xx, yy] = meshgrid (tx, ty);
 mesh (tx, ty, w');
 xlabel ('x');
 ylabel ('y');
 zlabel ('div(x,y)');
 title ('div test');

