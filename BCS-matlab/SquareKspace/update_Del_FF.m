function DiffDel=update_Del_FF()
global Del TT DD HBdG Nx Ny mu t hh kx ky Temp U0 DelHelicity
global Density Pol Ek UUk alpha ED UU QQ DensityHelicity DDH UTT UTTb 

Delold=Del;

Del0=sum(sum(Del))/(Nx*Ny);
DD=[0,Del0;-Del0 0];   


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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Helicity Basis
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [UTT ETT]=eig(TT);
    [ETT, order]=sort(diag(ETT),'descend');
    UTT=UTT(:,order);
    
%      v1=exp(-1i*angle(UTT(1,1)));
%      UTT(:,1)=v1*UTT(:,1);
%      
%      v1=exp(-1i*angle(UTT(2,2)));
%      UTT(:,2)=v1*UTT(:,2);
    
    [UTTb ETTb]=eig(TTb);
    
    UTTb=UTT;
    UTTb(1,2)=-UTTb(1,2);
    UTTb(2,1)=-UTTb(2,1);
        
    TT=UTT'*TT*UTT;
      
    DDH=UTT'*DD*conj(UTTb);
    
    TTb=UTTb'*TTb*UTTb;
        
    HBdG=(HBdG+HBdG')/4;
    
    HBdG=[ TT DDH ; DDH' -TTb.'];
    
    [UU,ED]=eig(HBdG);
    [ED, order]=sort(diag(ED),'descend');
    UU=UU(:,order);
    FEV=FD(ED,Temp);
    
    DelHelicity(i,j)=U0*UU(2,:)*(UU(4,:)'.*FEV(:));
    %DelHelicity(i,j)=U0*UU(1,:)*(UU(3,:)');
    
    DensityHelicity(1,i,j)=UU(1,:)*(UU(1,:)'.*FEV(:));
    DensityHelicity(2,i,j)=UU(2,:)*(UU(2,:)'.*FEV(:));   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Helicity Basis
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
end


DiffDel=sum(sum(abs(Del-Delold))/(Nx*Ny));





