function msweep(istart,iend)
global SopL SopR N Spsi Sops ti Nph mBt tisweep dt

x=3; % index of the site from left side for the correlation function (be careful x>=ti)
r=1; % distance between sites whose correlations in question


corr=zeros(1, N-2*ti-1);      %initial array for correlation function

for i=istart:-1:iend

 % Sx_i*Sx_{i+1} correlation
 Sxx=kron( kron(SopL{i,1},Sops{5}) , kron(Sops{5},SopR{N-i-2,1}) );
 Sxy=kron(       kron(SopL{i,1},Sops{5}) ,   -0.5i*kron(Sops{2}-Sops{3},SopR{N-i-2,1}) );
 Syx=kron( -0.5i*kron(SopL{i,1},Sops{2}-Sops{3}) , kron(Sops{5},        SopR{N-i-2,1}) );
 
 size(Sxy)
 size(Syx)

 psi=Spsi{i}; 

 %corr(N-i-ti-1)=psi'*(Sxy-Syx)*psi/4.;
 corr(N-i-ti-1)=psi(:,1)'*Sxx*psi(:,1)/4.;
 
 mx=psi'*( ...
           kron( kron(SopL{i,5},Sops{1}) , kron(Sops{1},SopR{N-i-2,1}) ) ...
          +kron( kron(SopL{i,1},Sops{5}) , kron(Sops{1},SopR{N-i-2,1}) ) ...
          +kron( kron(SopL{i,1},Sops{1}) , kron(Sops{5},SopR{N-i-2,1}) ) ...
          +kron( kron(SopL{i,1},Sops{1}) , kron(Sops{1},SopR{N-i-2,5}) ) ...
         )*psi/N;
 
end

 sum(corr)/(istart-iend+1)

dlmwrite(strcat('data/Heisenberg_N',int2str(N),'_mBt',int2str(mBt),'.dat'), [(istart:-1:iend)'  tisweep*ones(size(corr',1),1) corr'], 'delimiter', '\t', '-append');

disp('wrote data')

%file_id = fopen(strcat('data/gs_occupdist_N',int2str(N),'_Nph',int2str(Nph),'_mBt',int2str(mBt),'.dat'), 'a');
%fprintf(file_id, ' %i \n', ti);
%printf(' %i \n', ti)

