function energy_FF()

global Del TT Nx Ny mu t hh kk Temp U0 
global  Ek  band1 QQ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Energy per site
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Etot=0;
TT=0;
band1=zeros(Nx,Ny);
band2=zeros(Nx,Ny);
band3=zeros(Nx,Ny);
band4=zeros(Nx,Ny);
for i=1:Nx    
    for j=1:Ny
        kpQ(1)=kk(i,j,1)+QQ(1);     
        kpQ(2)=kk(i,j,2)+QQ(2);
        
         TT(1,1)=-mu-2*t*(cos(kpQ(1))+cos(0.5*( kpQ(1)+sqrt(3)*kpQ(2)))...
                                     +cos(0.5*(-kpQ(1)+sqrt(3)*kpQ(2))) )-hh(3);
         TT(2,2)=-mu-2*t*(cos(kpQ(1))+cos(0.5*( kpQ(1)+sqrt(3)*kpQ(2)))...
                                     +cos(0.5*(-kpQ(1)+sqrt(3)*kpQ(2))) )+hh(3);
       
    
        Etot=Etot+trace(TT)/2;    
        Etot=Etot+FD(Ek(:,i,j)',Temp)*Ek(:,i,j);   
        
        band1(i,j)=Ek(1,i,j);         
        band2(i,j)=Ek(2,i,j);        
        band3(i,j)=Ek(3,i,j);        
        band4(i,j)=Ek(4,i,j);
      
    end
end

Del0=sum(sum(Del))/(Nx*Ny); 
Etot=Etot/(Nx*Ny) +abs(Del0)^2/U0

%return

% hold on
% mesh(band1)
% mesh(band2)
% mesh(band3)
% mesh(band4)
% hold off











