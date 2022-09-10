function plotpol3()


global Nx Ny kx ky Pol


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
[U, V, W]=surfnorm(X,Y,polx);
surf(X,Y,polx)
hold on
quiver3(X,Y,polx,U,V,W)
colormap jet
colorbar
title('Polarization in X direction')
hold off

%%%%%%%%%%%%%%%
figure(2)
[U, V, W]=surfnorm(X,Y,poly);
surf(X,Y,poly)
hold on
quiver3(X,Y,poly,U,V,W)
colormap hot
colorbar
title('Polarization in Y direction')
hold off
%%%%%%%%%%%%%%%
figure(3)
polz=real(polz);
[U, V, W]=surfnorm(X,Y,polz);
surf(X,Y,polz)
hold on
quiver3(X,Y,polz,U,V,W)
colormap cool
colorbar
title('Polarization in Z direction')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%