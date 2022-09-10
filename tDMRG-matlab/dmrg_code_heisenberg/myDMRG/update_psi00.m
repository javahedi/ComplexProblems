function psi0=update_psi00(i,psi,block) 
global SOL SOR N mBt ms Id Spsi mBt ti

%mB=ms*mBt;
%mBr=ms*mBt;
%n01=norm(psi)
%sqrt(psi'*psi)
SOL{ti}=SOL{ti+1};
SOR{ti}=SOR{ti+1};
if block==-1
  OL=SOL{i+1};
  OR=SOR{N-i-2};
  mB=size(OL,1);
  mBr=size(OR,2)*ms;
  Mpsi=reshape(psi,mBr,mB);
  Mpsi=Mpsi.';
  Mpsi0=OL'*Mpsi*kron(Id,OR).';
  mB=size(OL,2);
  mBr=size(OR,1)*ms;

elseif block==+1

  OL=SOL{i};
  OR=SOR{N-i-1};
  mB=size(OL,2)*ms;
  mBr=size(OR,1);
  Mpsi=reshape(psi,mBr,mB);
  Mpsi=Mpsi.';
  Mpsi0=kron(OL,Id)*Mpsi*conj(OR);
  mB=size(OL,1)*ms;
  mBr=size(OR,2);
end

Mpsi0=Mpsi0.';

psi0=reshape(Mpsi0,mB*mBr,1);

%n02=norm(psi0)
%sqrt(psi0'*psi0)

