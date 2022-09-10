function DiffDel=update_Del_FF()

global Del TT DD HBdG Nx Ny mu t hh kx ky Temp U0
global Density Pol Ek UUk alpha ED UU QQ

Delold=Del;

Del0=sum(sum(Del))/(Nx*Ny);
DD=[0,Del0;-Del0 0];   

%v1=zeros(1,4);

for i=1:Nx
    
    kpQ(1)=kx(i)+QQ(1);
    %kpQ(1)=kpQ(1)-sign(kpQ(1))*floor( abs(kpQ(1))/pi )*2*pi;
    
    for j=1:Ny
    
    kpQ(2)=ky(j)+QQ(2);
    %kpQ(2)=kpQ(2)-sign(kpQ(2))*floor( abs(kpQ(2))/pi )*2*pi;
        
    TT(1,1)=-mu-2*t*(cos(kpQ(1))+cos(kpQ(2)))-hh(3);
    TT(1,2)=-hh(1)+1i*hh(2)+2*alpha*(-1i*sin(kpQ(1))+sin(kpQ(2)));
    TT(2,1)=-hh(1)-1i*hh(2)+2*alpha*(1i*sin(kpQ(1))+sin(kpQ(2)));
    TT(2,2)=-mu-2*t*(cos(kpQ(1))+cos(kpQ(2)))+hh(3);

    TTb(1,1)=-mu-2*t*(cos(-kx(i))+cos(-ky(j)))-hh(3);
    TTb(1,2)=-hh(1)+1i*hh(2)+2*alpha*(-1i*sin(-kx(i))+sin(-ky(j)));
    TTb(2,1)=-hh(1)-1i*hh(2)+2*alpha*(1i*sin(-kx(i))+sin(-ky(j)));
    TTb(2,2)=-mu-2*t*(cos(-kx(i))+cos(-ky(j)))+hh(3);
    

    HBdG=[ TT DD ; DD' -TTb.'];
  
    
    HBdG=(HBdG+HBdG')/4;

    [UU,ED]=eig(HBdG);
    [ED, order]=sort(diag(ED),'descend');
    UU=UU(:,order);
    FEV=FD(ED,Temp);

    %v1=diag(exp(-1i*angle(UU(1,:))));
    %UU=v1*UU;

    UUk(:,:,i,j)=UU;
    Ek(:,i,j)=ED;
    
            
    Del(i,j)=U0*UU(1,:)*(UU(4,:)'.*FEV(:));
    %Del2(i,j)=U0*UU(:,4)'*(FEV(:).*UU(:,1));
    %Del3(i,j)=U0*UU(2,:)*(UU(3,:)'.*FEV(:));

    Density(1,i,j)=UU(1,:)*(UU(1,:)'.*FEV(:));
    Density(2,i,j)=UU(2,:)*(UU(2,:)'.*FEV(:));
    
    Pol(1,i,j)=2*real( UU(2,:)*(UU(1,:)'.*FEV(:))  );
    Pol(2,i,j)=2*imag( UU(2,:)*(UU(1,:)'.*FEV(:))  );
    Pol(3,i,j)=Density(1,i,j)-Density(2,i,j);
    
    %topology(i,j) ; 
    
    end
end


DiffDel=sum(sum(abs(Del-Delold))/(Nx*Ny));





