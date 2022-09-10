function [psi, Es]=super_block(i,opLe,opRe,block)
global Jx Jz  Spsi  Nop  Nmid N
   
   Hs=kron(opLe{Nop},opRe{1})+kron(opLe{1},opRe{Nop}) ...
     +Jx*0.5*(kron(opLe{2},opRe{3})+kron(opLe{3},opRe{2}))...
     +Jz*kron(opLe{4},opRe{4});
  
   Hs=(Hs+Hs')/2;  % to make sure Hs is numerically Hermitian
  

   if block==0 
    [psi, Es]=eigs(Hs,1,'sa');
    Spsi{i+1}=psi;
   else
   % opts.v0=Spsi{i};
    opts.maxit=10000;    % maximum number of iterations
   % opts.disp=1;
    [psi2, Es2]= eigs(Hs,4,'sa',opts);
    %[psiv, Ev]=eig(Hs);
    %psi=psiv(:,1);
    %Es=Ev(1,1);
    psi=psi2(:,4);
    Es = Es2(4,4); 
   % Es2
    -(Es2(4,4)-Es2(3,3))/N
   end
   Es/N

 
