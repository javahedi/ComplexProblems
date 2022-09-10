function prepare_tsweep

global U J V0 ti  Nop
global Sops SopL SopR 

  [vl, vr]=potential(ti-1,1,0);
  Sops{Nop}=U*Sops{6};
  if ti==1
    SopL{ti,4}=vl*Sops{4};
    SopR{ti,4}=vr*Sops{4};
    SopL{ti,Nop}=U*SopL{ti,6}+V0*SopL{ti,4};             % block has 1 site, no hopping inside the block!
    SopR{ti,Nop}=U*SopR{ti,6}+V0*SopR{ti,4};
  else
    SopL{ti,4}=vl*kron(SopL{ti-1,1},Sops{4});
    SopR{ti,4}=vr*kron(Sops{4}, Sops{ti-1,1});
    SopL{ti,Nop}=-J*(kron(SopL{ti-1,3},Sops{2})+kron(SopL{ti-1,2},Sops{3}))+U*SopL{ti,6}+V0*SopL{ti,4};
    SopR{ti,Nop}=-J*(kron(Sops{2}, SopR{ti-1,3})+kron(Sops{3},SopR{ti-1,2}))+U*SopR{ti,6}+V0*SopR{ti,4};
end
