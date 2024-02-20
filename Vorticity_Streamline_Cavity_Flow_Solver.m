%%  VORTICITY/STREAM FUNCTION LID DRIVEN CAVITY FLOW SOLVER   JOE MOLVAR
clear;  close all
%%% GIVENS
Nx = 32;  L = 1; Wall_Velocity = 1;       % Nodes X; Domain Size; Velocity
rho = 1;   mu = 0.01;                     % Density; Dynamic Viscosity;
dt = 0.001;  maxIt = 50000;  maxe = 1e-7; % Time Step; Max iter; Max error
%%% SETUP 1D GRID
Ny = Nx;  h=L/(Nx-1);  x = 0:h:L;  y = 0:h:L;  
im = 1:Nx-2; i = 2:Nx-1; ip = 3:Nx;  jm = 1:Ny-2; j = 2:Ny-1; jp = 3:Ny;
%%% PRELOCATE MATRIXES
Vo = zeros(Nx,Ny);  St = Vo;  Vop = Vo;  u = Vo;  v = Vo;           
%%% SOLVE LOOP SIMILAR TO GAUSS-SIEDEL METHOD
for iter = 1:maxIt     
    %%% CREATE BOUNDARY CONDITIONS
    Vo(1:Nx,Ny) = -2*St(1:Nx,Ny-1)/(h^2) - Wall_Velocity*2/h; % Top
    Vo(1:Nx,1)  = -2*St(1:Nx,2)   /(h^2);                     % Bottom
    Vo(1,1:Ny)  = -2*St(2,1:Ny)   /(h^2);                     % Left
    Vo(Nx,1:Ny) = -2*St(Nx-1,1:Ny)/(h^2);                     % Right
    %%% PARTIALLY SOLVE VORTICITY TRANSPORT EQUATION
    Vop = Vo;
    Vo(i,j) = Vop(i,j) + ...
        (-1*(St(i,jp)-St(i,jm))/(2*h) .* (Vop(ip,i)-Vop(im,j))/(2*h)+...
            (St(ip,j)-St(im,j))/(2*h) .* (Vop(i,jp)-Vop(i,jm))/(2*h)+...
     mu/rho*(Vop(ip,j)+Vop(im,j)-4*Vop(i,j)+Vop(i,jp)+Vop(i,jm))/(h^2))*dt;
    %%% PARTIALLY SOLVE ELLIPTICAL VORTICITY EQUATION FOR STREAM FUNCTION
    St(i,j) = (Vo(i,j)*h^2 + St(ip,j) + St(i,jp) + St(i,jm) + St(im,j))/4;
    %%% CHECK FOR CONVERGENCE
    if iter > 10
        error = max(max(Vo - Vop))
        if error < maxe
            break;
        end
    end
 end
%%% CREATE VELOCITY FROM STREAM FUNCTION
u(2:Nx-1,Ny) = Wall_Velocity;
u(i,j) = (St(i,jp)-St(i,jm))/(2*h);  v(i,j) = (-St(ip,j)+St(im,j))/(2*h);
%%% PLOTS
cm = hsv(ceil(100/0.7));  cm = flipud(cm(1:100,:));
figure(1);  contourf(x,y,u',23,'LineColor','none');
title('U-velocity');  xlabel('x-location');  ylabel('y-location')  
axis('equal',[0 L 0 L]);  colormap(cm);  colorbar('westoutside');   
figure(2);  plot(y,u(round(Ny/2),:));
title('Centerline x-direction velocity');
xlabel('y/L');  ylabel('u/U');  axis('square');  xlim([0 L]);  grid on  
N = 1000;  xstart = max(x)*rand(N,1);  ystart = max(y)*rand(N,1);
[X,Y] = meshgrid(x,y);
figure(3);  h=streamline(X,Y,u',v',xstart,ystart,[0.1, 200]);
title('Stream Function');  xlabel('x-location');  ylabel('y-location')
axis('equal',[0 L 0 L]);  set(h,'color','k')