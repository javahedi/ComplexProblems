function init_mat()

global ms Sops N Id
global Nop SopL SopR J1 D1 HH DM

Sops=cell(6);
SopL=cell(N,Nop);
SopR=cell(N,Nop);

%initialization of matrices 

ms=2; % hilbert space dimension of a single site

Hs=zeros(ms,ms);%  H for a site
Id=eye(ms,ms);

sx=[0 1; 1 0];
sy=[0 -1i; 1i 0];
sz=[1 0; 0 -1];
sp=sx+1i*sy;
sm=sx-1i*sy;

%block operator initializations 
Sops{1}=Id;
Sops{2}=sp;
Sops{3}=sm;
Sops{4}=sz;
Sops{5}=sx;    
Sops{6}=Hs;

%block operator initializations for storing (left and right) 
SopL{2,1}=kron(Id,Id);
% spin ops of the right-most site of the left block
SopL{2,2}=kron(Id,sp);
SopL{2,3}=kron(Id,sm);
SopL{2,4}=kron(Id,sz);
% spin ops of the second-right site of the left block
SopL{2,5}=kron(sp,Id);
SopL{2,6}=kron(sm,Id);
SopL{2,7}=kron(sz,Id);
% total Sx
SopL{2,8}=kron(sx,Id)+kron(Id,sx);
% Hamiltonian
SopL{2,9}=J1*( 0.5*(kron(sp,sm)+kron(sm,sp))+D1*kron(sz,sz) ) +...
          DM*( 1i*0.5*(kron(sp,sm)-kron(sm,sp))) +...
             HH(1)*SopL{2,8};

SopR{2,1}=kron(Id,Id);
SopR{2,2}=kron(sp,Id);
SopR{2,3}=kron(sm,Id);
SopR{2,4}=kron(sz,Id);
SopR{2,5}=kron(Id,sp);
SopR{2,6}=kron(Id,sm);
SopR{2,7}=kron(Id,sz);
SopR{2,8}=kron(sx,Id)+kron(Id,sx);
SopR{2,9}=J1* ( 0.5*(kron(sp,sm)+kron(sm,sp))+ D1*kron(sz,sz) ) +...
          DM* ( 1i*0.5*(kron(sp,sm)-kron(sm,sp))) +...
           HH(1)*SopR{2,8};

