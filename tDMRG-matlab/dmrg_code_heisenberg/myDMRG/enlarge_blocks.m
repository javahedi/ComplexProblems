function [opLe, opRe]=enlarge_blocks(i,block)


global Sops Jx Jz hx  SopL SopR N Nop tisweep


if block==0
 opL=SopL(i,:);
 opR=SopR(i,:);
else
 opL=SopL(i,:);
 opR=SopR(N-i-2,:);
end

for i=1:Nop
  opLe{i}=kron(opL{1},Sops{i});
end
  
opLe{Nop-1}=opLe{Nop-1}+kron(opL{Nop-1},Sops{1});
opLe{Nop}=opLe{Nop}+kron(opL{Nop},Sops{1})+Jx*0.5*(kron(opL{2},Sops{3})+kron(opL{3},Sops{2}))+Jz*kron(opL{4},Sops{4})+hx*opLe{Nop-1};


for i=1:Nop
  opRe{i}=kron(Sops{i},opR{1});
end

opRe{Nop-1}=opRe{Nop-1}+kron(Sops{1},opR{Nop-1});
opRe{Nop}=opRe{Nop}+kron(Sops{1},opR{Nop})+Jx*0.5*(kron(Sops{2},opR{3})+kron(Sops{3},opR{2}))+Jz*kron(Sops{4},opR{4})+hx*opRe{Nop-1};


