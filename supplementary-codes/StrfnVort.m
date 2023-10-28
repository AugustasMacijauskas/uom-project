
%% Below is a Matlab code to integrate in time the vorticity equation


%   First change to $mu=0.1$, $nx=11$, $dt=0.035$, $dpic=0.025$ and $tfinal=0.525$.
%     The plot of the streamfunction should show the numerical instability.
%     Decreasing $dt$ to 0.3 nearly stabilises the result at this $tfinal$, while $dt=0.25$ works to $tfinal=3.0$.
%     This value of $dt$ is the marginal value, but that is for large grids, and small grids are slightly more stable.
%      Best to work at $dt=0.2*h*h/mu$

%  The code first calculates the boundary conditions to first-order, and then extrapolates these results to apply the conditions to second-order.
%    Use this second-order code.
%    Change to  $mu=0.1$, $nx=11$, $dpic=0.1$ and $tfinal=1.0$ and $dt=0.02$.
%    Find the value of $\omega(x=0.5,y=0.5,t=1.0)$.
%    Now decrease $dt$ to 0.01, 0.005, 0.0025 and 0.001.
%    Then with $nx=15$ try $dt$ = 0.01 (largest stable value), 0.005, 0.0025 and 0.001.
%    Finally with $nx=21$ try $dt=$ 0.005 (largest stable value), 0.0025 and 0.001.
%    Plot these results for $\psi(x=0.5,y=0.5,t=1.0)$ as a function of $dt$.
%    You could experiment by deleting the lines of the code that apply the boundary conditions at second-order.

%  Now set $tfinal=3.0$, $dpic=0.1$ and $dt=0.2*h*h/mu$ and obtain $\omega(x=0.5,y=0.5,t)$ to find how long it takes the vorticity to attain a steady value within 4 significant figures.
%     Compare your plots for the steady state of the streamfunction and the vorticity.

%  Gather results for different spatial resolutions $nx$ for the steady horizontal velocity $u$ at the mid-section $x=0.5 $, and plot on top of one another.

%  The code calculates the force on the top plate to second-order accuracy. Find the steady force for different spatial resolutions, $nx=11$, 15, 21, 29 and 41. Plot the force as a function of the grid size $h$.

% Change the top slip boundary condition from $u=\sin^2\pi x$ to $u=1$, i.e change
%%         w(i,nx)= -(sin(pi*(i-1)*h)*sin(pi*(i-1)*h)... 
%%   to    w(i,nx) = -1...
%  Find the  force on the top plate for various spatial resolutions, say $nx= 10$, 14, 20, 28 and 40. Show that the force diverges as the resolution increases as $F = 4.32\ln(1/h) - 3.75$.

 %% Streamfunction vorticity code
 clear all %clear all variables
 clc % clear output screen

 %Constants
 dpic	  = 0.5;
 tfinal = 3.0;
 mu	  = 0.1;  %i.e. Re=10
 nx	  = 10+1;
 ny	  = nx;  %different nx & ny gives non-square
 itmax    = 2*nx;
   
 w = zeros(nx,ny); % form a matrix of zeros (1..nx, 1..ny) 
 p = w; % w = vorticty, p = streamfunction
 wt = w; % wt = dw/dt
   
 %Initialization
 h = 1/(nx-1); %h = dx (=dy); so x-size is 1, y-size is ny/nx
 t = 0;
 r = 2/(1+ pi/(nx-1)); %SOR parameter
 dt = 0.2*h*h/mu; 
   
 for i= 1:nx
   for j=1:ny
     p(i,j) = 0;
     w(i,j) = 0;
   end
 end
   
 %% Main program starts here
 tpic = t + dpic - 0.5*t;
 while(t<tfinal)
   while(t<tpic)

     % Main Poisson solver: solve nabla^2 p = - w
     for it = 1:itmax % iterations start here
       for i = 2:nx-1 
         for j = 2:ny-1
           p(i,j) = (1-r)*p(i,j) + r*0.25*(p(i-1,j) + p(i+1,j) + ...
                  p(i,j-1) + p(i,j+1) + h*h*w(i,j) );  
         end
       end
     end % iterations end here
       
     % Time-step vorticity
           
     % BC 1st order
     for i= 2:nx-1
       w(i,1) =  -((p(i,2)-p(i,1))/h - 0)/(0.5*h);
       w(i,ny)= -(sin(pi*(i-1)*h)*sin(pi*(i-1)*h)... 
                    -(p(i,ny)-p(i,ny-1))/h)/(0.5*h);
     end
     for j= 2:ny-1
       w(1,j) = ((-(p(2,j)-p(1,j))/h) - 0)/(0.5*h);
       w(nx,j) = (0 - (-(p(nx,j)-p(nx-1,j))/h))/(0.5*h);
     end
           
     %BC 2nd order correction
     for i = 2:nx-1
       w(i,1) =  (4*w(i,1)  - w(i,2))/3;
       w(i,ny) = (4*w(i,ny) - w(i,ny-1))/3;
     end
     for j = 2:ny-1
       w(1,j) =  (4*w(1,j)  - w(2,j))/3;
       w(nx,j) = (4*w(nx,j) - w(nx-1,j))/3;
     end
     % end of bcs
        
     for i = 2:nx-1
       for j = 2:ny-1  
         wt(i,j) =-(p(i,j+1)-p(i,j-1))*(w(i+1,j)-w(i-1,j))/(4*h*h)...
            + (p(i+1,j)-p(i-1,j))*(w(i,j+1)-w(i,j-1))/(4*h*h)...
            + mu*(w(i+1,j)+w(i-1,j)+w(i,j+1)+w(i,j-1)-4*w(i,j))/(h*h);
       end
     end
        
     for i = 2:nx-1
       for j = 2:ny-1
         w(i,j) = w(i,j) + dt*wt(i,j);       
       end
     end
     t = t + dt;
        
   end %endwhile t<tpic

   formatSpec2='t=%10.5f; w(0.5, 0.5)=%10.5f\n ';
   fprintf(formatSpec2, t, w(round(nx/2),round(ny/2)));
   tpic = t + dpic - 0.5*dt;
 end %endwhile t<tfinal
    
 %Plot 'p', streamfunction 
 figure(1)   
 tx = linspace (0, (nx-1)*h, nx)';
 ty = linspace (0, (ny-1)*h, ny)';    
 [xx, yy] = meshgrid (tx, ty);
 mesh (tx, ty, p');
 xlabel ("x");
 ylabel ("y");
 zlabel ("p(x,y)");
 title ("Streamfunction Plot");
 % End of plot 'p'
    
 %Plot 'w' vorticity
 figure(2)    
 tx = linspace (0, (nx-1)*h, nx)';
 ty = linspace (0, (ny-1)*h, ny)';
 [xx, yy] = meshgrid (tx, ty);
 mesh (tx, ty, w');
 xlabel ("x");
 ylabel ("y");
 zlabel ("w(x,y)");
 title ("Vorticity Plot");
 % End of plot 'w'
        
 % midu: 'u' along x=0.5
 i = round(nx/2);
 for j = 2:ny
   formatSpec3='y=%10.5f; u=%10.5f\n ';
   fprintf(formatSpec3, h*(j-1-0.5), (p(i,j)-p(i,j-1))/h);
 end

 % Force calculation   
 q = 0;
 for i = 2:nx-1 
   q = q + (2*p(i,ny) - 5*p(i,ny-1) + 4*p(i,ny-2) -p(i,ny-3))/h;
   %1st order: (p[i,ny] - 2*p[i,ny-1] + p[i,ny-2])/h;
   %2nd order:(2*p[i,ny] - 5*p[i,ny-1] + 4*p[i,ny-2] -p[i,ny-3])/h;
 end
 formatSpec1='t=%10.4f; q=%10.5f\n;'; 
 fprintf(formatSpec1, t, q );
 % End of Force calculation
       
 % Final output
 formatSpec4='h=%10.5f; w(0.5,0.5)=%10.5f\n ';
 fprintf(formatSpec4, h, w(round(nx/2),round(ny/2)));
