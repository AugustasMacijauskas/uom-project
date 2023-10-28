
%Below is a Matlab code to solve the Poisson equation
%\[
% \nabla^2\psi = - 2\pi^2\sin\pi x\sin\pi y
%   \quad\text{in } 0\le x\le 1, 0 \le y \le 1,
%\]
%with boundary condition $\psi=0$ on $x=0$ and 1 and $y=0$ and 1.

%%   First, examine the graph of the solution for $\psi$ to see if is sensible.

%%   Second, see how the solution for the middle of the grid $\psi(0.5,0.5)$ converges with the number of iterations. 
%   Change the value of the SOR parameter $r$ from its preset optimal value.

%%   Find how everything depends on the number of points $nx$.

%%  The final test is to compare the numerical answer at different resolutions $nx$ with the exact analytic solution
%\[
% \psi = \sin\pi x\sin\pi y
%\]
%The code calculates the maximum error over the grid.
%    Examine how this error depends on $nx$.
%     Are the local errors $O(\Delta x^2)$.

 %% Poisson Test Program
  clear all  %clear all variables
  clc % clear output screen

  nx = 1+20; % Matlab's index for matrices starts from 1
  ny = nx; %different nx & ny gives non-square
  w	= zeros(nx,ny); % form a matrix of zeros (1..nx, 1..ny)
  p = w; % w = voriticty, p = streamfunction

  h = 1/(nx-1); %h = dx (=dy),
   % so x-size is 1, y-size is ny/nx

  % r = SOR parameter, below is optimal, try other values
  r = 2/(1+pi/(nx-1)); 

 %% Initiaze the values for p and w
  for i = 1:nx
    for j =1:ny
      p(i,j) = 0;
      w(i,j) = 2*pi*pi*sin(pi*(i-1)*h)*sin(pi*(j-1)*h);
    end
  end

 %%Iterations start here
  % max set to 4*nx. Try other values
  for it = 1:4*nx

    for i = 2:nx-1
      for j = 2:ny-1
        q = p(i-1,j) + p(i+1,j) + p(i,j-1) + p(i,j+1) -4*p(i,j) + h*h*w(i,j);
        p(i,j) =  p(i,j) +  r*0.25*q;
      end
    end 
     
    %Output during iterations
    formatSpec1='iteration=%d; psi(0.5,0.5) =%10.5f\n ';
    fprintf(formatSpec1, it, p(floor(nx/2),floor(nx/2))); 
    
  end %iterations

 %%er = error cf analytic answer for square box
  er = 0;
    for i = 1:nx
      for j = 1:ny
        if (abs(p(i,j)-sin(pi*(i-1)*h)*sin(pi*(j-1)*h)) > er) 
          er = abs(p(i,j)-sin(pi*(i-1)*h)*sin(pi*(j-1)*h));
        end
      end
    end

 %Output error
  formatSpec1='iteration=%d; er=%12.7f; p(0.5,0.5) =%10.5f\n ';
  fprintf(formatSpec1, it, er, p(floor(nx/2),floor(nx/2)));

 %%Plot (3D) solution
  tx = ty = linspace (0, (nx-1)*h, nx)';
  [xx, yy] = meshgrid (tx, ty);
  mesh (tx, ty, p);
  xlabel ("x");
  ylabel ("y");
  zlabel ("phi(x,y)");
  title ("Poisson Plot");


