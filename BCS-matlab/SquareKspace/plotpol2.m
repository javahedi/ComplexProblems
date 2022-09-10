function plotpol2()


global Nx Ny kx ky Pol polx poly polz


polx=zeros(Nx,Ny);
poly=zeros(Nx,Ny);
polz=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        polx(i,j)=Pol(1,i,j);
        poly(i,j)=Pol(2,i,j);
        polz(i,j)=Pol(3,i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(ky,kx);

figure(1) 
%[U V]=surfnorm(polx);
%U=Pol(1,:,:);
U=polx;
%V=Pol(2,:,:);
V=poly;
%contour(X,Y,polx)
%hold on
quiver(X,Y,U,V)
colormap cool
colorbar
title('Polarization in X direction')
%hold off
axis([-pi pi -pi pi])
return
%%%%%%%%%%%%%%%
figure(2)
[U V]=surfnorm(poly);
%contour(X,Y,poly)
%hold on
quiver(X,Y,U,V,'color',[0 1 0])
colormap cool
colorbar
title('Polarization in y direction')
%hold off
axis([-pi pi -pi pi])
%%%%%%%%%%%%%%%
figure(3)
[U V]=surfnorm(polz);
%contour(X,Y,polz)
%hold on
quiver(X,Y,U,V)
colormap cool
title('Polarization in z direction')
colorbar
%hold off
axis([-pi pi -pi pi])
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example

% [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%  z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%  contour(x,y,z), 
%  hold on
%  indices = 11:21;
% quiver(x(:,indices),y(:,indices),px(:,indices),py(:,indices),'color',[0 0 1]);
% quiver(x(:,1:10),y(:,1:10),px(:,1:10),py(:,1:10),'color',[1 0 0]);





