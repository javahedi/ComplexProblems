function Etot=energy_FF()

global Del TT Nx Ny mu t hh kx ky Temp U0 
global  Ek  QQ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Energy per site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Etot=0;
TT=0;
for i=1:Nx    
    for j=1:Ny
        TT(1,1)=-mu-2*t*(cos(kx(i)+QQ(1))+cos(ky(j)+QQ(2)) )-hh(3);
        TT(2,2)=-mu-2*t*(cos(kx(i)+QQ(1))+cos(ky(j)+QQ(2)) )+hh(3);
    
        Etot=Etot+trace(TT)/2;    
        Etot=Etot+FD(Ek(:,i,j)',Temp)*Ek(:,i,j);  
        
     end
end

Del0=sum(sum(Del))/(Nx*Ny); 
Etot=Etot/(Nx*Ny) +abs(Del0)^2/U0;













