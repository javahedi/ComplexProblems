clear all
clc

% Set parameter
L  = 10;
J = 1;
deltaI = 0.5;
deltaF = 5;
hz = 0;


[Sx,Sy,Sz]=spin_spin(L);

% Step 1: creating initial hamiltonian
H = sparse(0);
for i=1:L
    j=mod(i,L)+1;
    H = H + J*(deltaI*cell2mat(Sx(i))*cell2mat(Sx(j))+...
               cell2mat(Sy(i))*cell2mat(Sy(j))+...
               cell2mat(Sz(i))*cell2mat(Sz(j))) + hz * cell2mat(Sz(i));

end


% Step 2: evaluate initial GS
%[v0,d0]=eigs(H,1,"sa");
[v0,d0]=eig(H);
clear H;
fprintf(' Ground state energy E0/N: %12.8f\n',d0(1)/L)


# Step 3: creating final hamiltonian
H = sparse(0);
for i=1:L
    j=mod(i,L)+1;
    H = H + J*(deltaF*cell2mat(Sx(i))*cell2mat(Sx(j))+...
               cell2mat(Sy(i))*cell2mat(Sy(j))+...
               cell2mat(Sz(i))*cell2mat(Sz(j))) + hz * cell2mat(Sz(i));
    
end

[v,d]=eig(H);
clear H;
cn = v0(:,1)'*v;
en=diag(d);
clear v; v0;

clear Sx; Sy; Sz;




file1 = fopen('LE_Fig1a.txt','w');
time = linspace(0,5,500);

for t= time
%    disp(t)
    LE = abs(sum(abs(cn).^2*exp(-1j*en*t))).^2;
    rate = -log(LE)/L;
    fprintf(file1,'%6.4f\t%12.8f\n',t,rate);
end
fclose(file1);

