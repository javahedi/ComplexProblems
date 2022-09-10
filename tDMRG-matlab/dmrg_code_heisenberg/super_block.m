function [rhoSB, psi, Es]=super_block(i,opLe,opRe,block)
global Jx Jz  Spsi  Nop N ti mB mBr ms mBt
   
   Hs=kron(opLe{Nop},opRe{1})+kron(opLe{1},opRe{Nop}) ...
     +Jx*0.5*(kron(opLe{2},opRe{3})+kron(opLe{3},opRe{2}))...
     +Jz*kron(opLe{4},opRe{4});
  
   Hs=(Hs+Hs')/2;  % to make sure Hs is numerically Hermitian

   if block==0 
    [psi, Es]=eigs(Hs,2,'sa');
    Spsi{i+1}=psi;
    lB=sqrt(size(psi,1));
    %rhoSB=reshape(psi,sqrt(size(psi,1)),sqrt(size(psi,1)));
    rhoSB=zeros(lB,lB);
    for k=1:2
      rhoSB=rhoSB+0.5*reshape(psi(:,k),lB,lB);
    end
   else
    opts.v0=Spsi{i};
    opts.v0=opts.v0(:,1);
    opts.maxit=1000;    % maximum number of iterations
   % opts.disp=1;
    [psi, Es]= eigs(Hs,2,'sa',opts);

    (Es(2,2)-Es(1,1))/N;

    % Targeting density matrix rhoSB
    mB=ms*mBt;
    mBr=ms*mBt;

    if block==-1
        if i==ti
            mB=ms*ms^ti;
            mBr=size(psi,1)/mB;
        elseif i==N-ti-2  
            mBr=ms*ms^ti;
            mB=size(psi,1)/mBr;
        end
    elseif block==+1
        if i==ti
            mBr=ms*ms^ti;
            mB=size(psi,1)/mBr;
        elseif i==N-ti-2  
            mB=ms*ms^ti;
            mBr=size(psi,1)/mB;
        end
    end
    
    
    rhoSB=zeros(mBr,mB);
    for k=1:2
      rhoSB=rhoSB+0.5*reshape(psi(:,k),mBr,mB);
    end
   end
   
   disp('Eg=====================')
   Es/N
   disp('Sx=====================')
   psi'*( kron(opLe{5},opRe{1}) + kron(opLe{1},opRe{5}) )*psi/N
end

 
