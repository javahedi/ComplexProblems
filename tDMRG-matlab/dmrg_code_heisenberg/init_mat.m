function init_mat()

global ms Sops N Id
global Nop mBt SopL SopR 


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
Sops{Nop}=Hs;
%block operator initializations for storing (left and right) 
SopL(1,:)=Sops;
SopR(1,:)=Sops;

