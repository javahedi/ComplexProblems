function [O rho]=den_mat(psi0,ib,block)  
global ms mBt SOL SOR ti  N Nmid isweep

mB=ms*mBt;
mBr=ms*mBt;

if block==-1
  if ib==ti
      mB=ms*ms^ti;
      mBr=size(psi0,1)/mB;
  elseif ib==N-ti-2  
      mBr=ms*ms^ti;
      mB=size(psi0,1)/mBr;
  end
elseif block==+1
  if ib==ti
      mBr=ms*ms^ti;
      mB=size(psi0,1)/mBr;
  elseif ib==N-ti-2  
      mB=ms*ms^ti;
      mBr=size(psi0,1)/mB;
  end
end
phi=reshape(psi0,mBr,mB);

if block==-1
  rhoL=phi.'*conj(phi);
  rhoL=(rhoL+rhoL')/2;
  [evec,eval]=eig(rhoL);
  [eval,order]=sort(diag(eval),'descend');
  evec=evec(:,order);
  O=evec(:,1:mBt);
  SOL{ib+1}=O;
  rho=rhoL;
  
%  if ib==Nmid-1
%    mBterror=sum((eval(mBt+1:end)));      
%    sumevalrhoL=sum(eval(1:mBt));                %sum of most weighted eigenvalues of rhoL 
%    data_eig=[eval sumevalrhoL*ones(size(eval,1),1)  mBterror*ones(size(eval,1),1) mBt*ones(size(eval,1),1)   isweep*ones(size(eval,1),1) ];
%    dlmwrite(strcat('data/rdm_eig_N_',int2str(N),'_mBt_',int2str(mBt),'.dat'), data_eig,'-append','roffset',1,'delimiter',' ', 'precision', 15);
%  end
  
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



