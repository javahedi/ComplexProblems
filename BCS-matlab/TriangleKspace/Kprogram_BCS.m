clc
clear all
format long
global Del TT DD HBdG Nx Ny N mu t hh Temp U0 
global Density Pol Ek UUk alpha GGk kk
global e1 e2 g1 g2

tic

Nx=8;                   % number of mesh over FBZ
Ny=32;
N=Nx*Ny;
e1=[1 0];             % lattice space set as constat a=1
e2=0.5*[1 sqrt(3)];
g1=(4*pi/(2*sqrt(3)*Nx))*[sqrt(3) -1];
g2=(4*pi/(  sqrt(3)*Ny))*[0 1];

if( mod(Nx,2)==1 )
    nx=(-(Nx-1)/2:1:Nx/2);
else
    nx=(-Nx/2+1:1:Nx/2);
end
if( mod(Ny,2)==1 )
    ny=(-(Ny-1)/2:1:Ny/2);
else
    ny=(-Ny/2+1:1:Ny/2);
end

kk=zeros(Nx,Ny,2);

for i1=1:Nx
    for i2=1:Ny
        kk(i1,i2,:)= nx(i1)*g1 + ny(i2)*g2;
    end
end



%kx=linspace(0,2*(Nx-1)*pi/Nx,Nx);   % mesh over kx in the FBZ
%kx=[ kx(ceil(Nx/2+.5)+1:end)-2*pi kx(1:ceil(Nx/2+.5)) ];
%ky=linspace(0,2*(Ny-1)*pi/Ny,Ny);   % mesh over ky in the FBZ
%ky=[ ky(ceil(Ny/2+.5)+1:end)-2*pi ky(1:ceil(Ny/2+.5)) ];

HBdG=zeros(4,4);        % Initial BdG hamiltonian
TT=zeros(2,2);          % Initial first quarter of HBdG matrix
DD=zeros(2,2);          % Initial delta
Ek=zeros(4,Nx,Ny);        % Energy spectrum 
UUk=zeros(4,4,Nx,Ny);     %
GGk=zeros(2,2,Nx,Ny);
mu=0;                   % Chemicsl potential
hh=zeros(1,3);          % Magnetic field  hx,hy,hz
alpha=1;                  % Spin-Orbit intersction
t=1;                    % hopping parameter
U0=-4;
Del=ones(Nx,Ny);          %Initial guess for gap order parameter
Density=zeros(2,Nx,Ny);
Pol=zeros(3,Nx,Ny);
hh(1:1)=0.5;

Temp=1.e-3;
DelEps=1.e-3;
itmax=100;

it=0;
DiffDel=2*DelEps;
while it<=itmax && DiffDel>DelEps
    it=it+1;
    DiffDel=update_Del_BCS;     
end
it



 %topology2_BCS();



toc

disp('Delta ')
sum(sum(Del))/(Nx*Ny)


% disp('Density ')
% sum(Density,2)/(Nx*Ny)
% disp('Ntot ')
% Ntot=sum(sum(Density))
% 

%disp('Polarization ')
%sum(sum(Pol,3),2)/(Nx*Ny)

 

%disp('Gamma ')
%Gamma=sum(sum(GGk,4),3)%*(2*pi)
%sum(Gamma)

%plot(kk,Ek,'.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Energy per site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energy_BCS()
 

 