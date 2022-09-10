function msweep(istart,iend)
global SopL SopR N Spsi Sops ti Nph mBt tisweep dt Nmid

x=3; % index of the site from left side for the correlation function (be careful x>=ti)
r=1; % distance between sites whose correlations in question


corr=zeros(1, N-2*ti-1);      %initial array for correlation function

for i=istart:-1:iend

  %if  x>=i || i>=x+r         % two observables are in the same block, they should be in the same basis
   
  %  opRe=kron(SopL{x+r,4},Sops{1});
  %  opLe=kron(SopL{x,4}, Sops{1});

  %elseif x<i && i<x+r        % two observables are in different blocks

  % opRe=kron(Sops{1},SopR{x+r,4});
  % opLe=kron(SopL{x,4},Sops{1});

  %end

 % Sx_i*Sx_{i+1} correlation
 Sxx=kron( kron(SopL{i,1},Sops{5}) , kron(Sops{5},SopR{i,1}) );
 
 psi=Spsi{i}; 

 %[O, rho]=den_mat(psi, i, -1);   % rhoL
 %corr(N-i-ti-1)=trace(rho*opLe*opRe);

 corr(N-i-ti-1)=psi'*Sxx*psi;
 
end

 %sum(corr)/(istart-iend+1);
 
mx=Spsi{Nmid}'*(kron(kron(SopL{Nmid-1,5}, Sops{1}), kron(Sops{1}, SopR{Nmid-1, 5})))*Spsi{Nmid}/N


dlmwrite(strcat('data/Heisenberg_N',int2str(N),'_mBt',int2str(mBt),'.dat'), [(istart:-1:iend)'  tisweep*ones(size(corr',1),1) corr'], 'delimiter', '\t', '-append');

%file_id = fopen(strcat('data/gs_occupdist_N',int2str(N),'_Nph',int2str(Nph),'_mBt',int2str(mBt),'.dat'), 'a');
%fprintf(file_id, ' %i \n', ti);
%printf(' %i \n', ti)

