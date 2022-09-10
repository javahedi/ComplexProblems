function [dHQx, dHQy]=dHQ()

global TT Nx Ny  t kx ky Temp 
global  Ek QQ 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Derivation respect to Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dHQx=0;

for i=1:Nx    
    for j=1:Ny
        TT(1,1)=2*t*sin(kx(i)+QQ(1));
        TT(2,2)=2*t*sin(kx(i)+QQ(1));
    
        dHQx=dHQx+trace(TT)/2;    
        dHQx=dHQx+FD(Ek(:,i,j)',Temp)*Ek(:,i,j);        
      
    end
end

dHQy=0;
for i=1:Nx    
    for j=1:Ny
        TT(1,1)=2*t*sin(ky(j)+QQ(2));
        TT(2,2)=2*t*sin(ky(j)+QQ(2));
    
        dHQy=dHQy+trace(TT)/2;    
        dHQy=dHQy+FD(Ek(:,i,j)',Temp)*Ek(:,i,j);         
      
    end
end
