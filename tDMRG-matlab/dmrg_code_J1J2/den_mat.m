function [O rho]=den_mat(rhoSB,ib,block)  
global ms mBt SOL SOR ti  N Nmid isweep

phi=rhoSB;

if block==-1
  rhoL=phi.'*conj(phi);
  rhoL=(rhoL+rhoL')/2;
  [evec,eval]=eig(rhoL);
  [eval,order]=sort(diag(eval),'descend');
  evec=evec(:,order);
  O=evec(:,1:mBt);
  SOL{ib+1}=O;
  rho=rhoL;
    
elseif block==1
  rhoR=phi*phi';
  rhoR=(rhoR+rhoR')/2;
  [evec,eval]=eig(rhoR);
  [eval,order]=sort(diag(eval),'descend');
  evec=evec(:,order);
  O=evec(:,1:mBt);
  SOR{ib+1}=O;
  rho=rhoR;
end



