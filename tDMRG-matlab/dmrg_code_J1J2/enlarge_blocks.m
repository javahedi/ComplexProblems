function [opLe, opRe]=enlarge_blocks(i,block)


global Sops J1 J2 D1 D2 HH  SopL SopR N Nop tisweep


if block==0
 opL=SopL(i,:);
 opR=SopR(i,:);
else
 opL=SopL(i,:);
 opR=SopR(N-i-2,:);
end

% add site and update identity and right-most of left block
for i=1:4
  opLe{i}=kron(opL{1},Sops{i});
end

% update second-right spin operators
for i=2:4
  opLe{i+3}=kron(opL{i},Sops{1});
end

opLe{8}=kron(opL{1},Sops{5}) + kron(opL{8},Sops{1});

opLe{9}=kron(opL{9},Sops{1}) ...
        + J1* ( 0.5*( kron(opL{2},Sops{3})+kron(opL{3},Sops{2}) )...
        + D1*kron(opL{4},Sops{4}) )...
        + HH(1)*(opLe{2}+opLe{3})/2 ...  % kron(opL{1},Sops{5})...
        + J2* ( 0.5*( kron(opL{5},Sops{3})+kron(opL{6},Sops{2}) )...
        + D2*kron(opL{7},Sops{4}) );

for i=1:4
  opRe{i}=kron(Sops{i},opR{1});
end

for i=2:4
  opRe{i+3}=kron(Sops{1},opR{i});
end

opRe{8}=kron(Sops{5},opR{1}) + kron(Sops{1},opR{8});

opRe{9}=kron(Sops{1},opL{9}) ...
        + J1* ( 0.5*( kron(Sops{3},opR{2})+kron(Sops{2},opR{3}) )...
        + D1*kron(Sops{4},opR{4}) )...
        + HH(1)*(opRe{2}+opRe{3})/2 ...
        + J2* ( 0.5*( kron(Sops{3},opR{5})+kron(Sops{2},opR{6}) )...
        + D1*kron(Sops{4},opR{7}) ) ;

