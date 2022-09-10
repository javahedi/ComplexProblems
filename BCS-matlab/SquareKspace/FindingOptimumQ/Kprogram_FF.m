function [dHQx, dHQy]=Kprogram_FF(QQi)


global Del TT DD HBdG Nx Ny N mu t hh kx ky Temp U0 
global Density Pol Ek UUk alpha GGk QQ 

tic

Nx=16;                   % number of mesh over FBZ
Ny=16;
N=Nx*Ny;
kx=linspace(0,2*(Nx-1)*pi/Nx,Nx);   % mesh over kx in the FBZ
kx=[ kx(ceil(Nx/2+.5)+1:end)-2*pi kx(1:ceil(Nx/2+.5)) ];
ky=linspace(0,2*(Ny-1)*pi/Ny,Ny);   % mesh over ky in the FBZ
ky=[ ky(ceil(Ny/2+.5)+1:end)-2*pi ky(1:ceil(Ny/2+.5)) ];
HBdG=zeros(4,4);        % Initial BdG hamiltonian
TT=zeros(2,2);          % Initial first quarter of HBdG matrix
DD=zeros(2,2);          % Initial delta
Ek=zeros(4,Nx,Ny);      % Energy spectrum 
UUk=zeros(4,4,Nx,Ny);   % EigenVector
GGk=zeros(2,2,Nx,Ny);   % Vector for Topological Number
mu=0;                   % Chemicsl potential
hh=zeros(1,3);          % Magnetic field  hx,hy,hz
alpha=0.75;             % Spin-Orbit intersction
t=1;                    % hopping parameter
U0=-4;
Del=ones(Nx,Ny);        %Initial guess for gap order parameter
Density=zeros(2,Nx,Ny);
Pol=zeros(3,Nx,Ny);
hh(1:2)=0.0;


QQ=QQi;

Temp=1.e-3;
DelEps=1.e-3;
itmax=100;

it=0;
DiffDel=2*DelEps;
while it<=itmax && DiffDel>DelEps
    it=it+1
    DiffDel=update_Del_FF;     
end
it

[dHQx, dHQy]=dHQ();

%topology2_FF();

toc

disp('Delta ');
sum(sum(Del))/(Nx*Ny);

% disp('Density ')
% sum(Density,2)/(Nx*Ny)
% disp('Ntot ')
Ntot=sum(sum(sum(Density,3),2));
% 
disp('Polarization ');
sum(sum(Pol,3),2)/(Nx*Ny);
% 

disp('Gamma ');
Gamma=sum(sum(GGk,4),3);%*(2*pi)
sum(Gamma);

%plot(kk,Ek,'.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Energy per site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%energy_FF()
 