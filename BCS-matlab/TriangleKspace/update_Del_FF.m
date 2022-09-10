function DiffDel=update_Del_FF()
global Del TT DD HBdG Nx Ny mu t hh kk Temp U0
global Density Pol Ek UUk alpha ED UU QQ DelHelicity

Delold=Del;

Del0=sum(sum(Del))/(Nx*Ny);
DD=[0,Del0;-Del0 0];   


for i=1:Nx      
    for j=1:Ny
        kpQ(1)=kk(i,j,1)+QQ(1);     
        kpQ(2)=kk(i,j,2)+QQ(2);   
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
                                  
    TT(1,1)=-mu-2*t*(cos(kpQ(1))+cos(0.5*( kpQ(1)+sqrt(3)*kpQ(2)))...
                                +cos(0.5*(-kpQ(1)+sqrt(3)*kpQ(2))) )-hh(3);
    
                         
    TT(1,2)=-hh(1)+1i*hh(2)-1i*2*alpha*(sin(kpQ(1))+0.5*sin(0.5*( kpQ(1)+sqrt(3)*kpQ(2)))...
                                                   -1i*0.5*sqrt(3)*sin(0.5*(kpQ(1)+sqrt(3)*kpQ(2)))...
                                                   -0.5*sin(0.5*( -kpQ(1)+sqrt(3)*kpQ(2)))...
                                                   -1i*0.5*sqrt(3)*sin(0.5*(-kpQ(1)+sqrt(3)*kpQ(2))));
    %TT(2,1)=conj(TT(1,2)) ;              
    TT(2,1)=-hh(1)-1i*hh(2)+1i*2*alpha*(sin(kpQ(1))+0.5*sin(0.5*( kpQ(1)+sqrt(3)*kpQ(2)))...
                                                  +1i*0.5*sqrt(3)*sin(0.5*(kpQ(1)+sqrt(3)*kpQ(2)))...
                                                  -0.5*sin(0.5*( -kpQ(1)+sqrt(3)*kpQ(2)))...
                                                  +1i*0.5*sqrt(3)*sin(0.5*(-kpQ(1)+sqrt(3)*kpQ(2))));
     %TT(2,2)=TT(1,1)+2*hh(3);                  
    TT(2,2)=-mu-2*t*(cos(kpQ(1))+cos(0.5*( kpQ(1)+sqrt(3)*kpQ(2)))...
                                 +cos(0.5*(-kpQ(1)+sqrt(3)*kpQ(2))) )+hh(3);
                               
                               
                               
                               
                               
      kpQ(1)=kk(i,j,1);     
      kpQ(2)=kk(i,j,2);           
                           
     TTb(1,1)=-mu-2*t*(cos(-kpQ(1))+cos(0.5*(-kpQ(1)-sqrt(3)*kpQ(2)))...
                                 +cos(0.5*( kpQ(1)-sqrt(3)*kpQ(2))) )-hh(3);
                           
                   
     TTb(1,2)=-hh(1)+1i*hh(2)-1i*2*alpha*(sin(-kpQ(1))+0.5*sin(0.5*( -kpQ(1)-sqrt(3)*kpQ(2)))...
                                                      -1i*0.5*sqrt(3)*sin(0.5*(-kpQ(1)-sqrt(3)*kpQ(2)))...
                                                      -0.5*sin(0.5*( kpQ(1)-sqrt(3)*kpQ(2)))...
                                                      -1i*0.5*sqrt(3)*sin(0.5*(kpQ(1)-sqrt(3)*kpQ(2))));
                      
    TTb(2,1)=-hh(1)-1i*hh(2)+1i*2*alpha*(sin(-kpQ(1))+0.5*sin(0.5*( -kpQ(1)-sqrt(3)*kpQ(2)))...
                                                     +1i*0.5*sqrt(3)*sin(0.5*(-kpQ(1)-sqrt(3)*kpQ(2)))...
                                                     -0.5*sin(0.5*( kpQ(1)-sqrt(3)*kpQ(2)))...
                                                     +1i*0.5*sqrt(3)*sin(0.5*(kpQ(1)-sqrt(3)*kpQ(2))));
                                                  
    TTb(2,2)=-mu-2*t*(cos(-kpQ(1))+cos(0.5*(-kpQ(1)-sqrt(3)*kpQ(2)))...
                                  +cos(0.5*( kpQ(1)-sqrt(3)*kpQ(2))) )+hh(3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    

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
    
    [UTT, ETT]=eig(TT);
    [UTTb, ETTb]=eig(TTb);
    
    UT=[UTT zeros(2); zeros(2) UTTb'];
            
    
    
    
    Del(i,j)=U0*UU(1,:)*(UU(4,:)'.*FEV(:));
    DelHelicity(i,j)=(UT(1,:)*UU(:,:))*((UT(4,:)*UU(:,:))'.*FEV(:));
    %Del2(i,j)=U0*UU(:,4)'*(FEV(:).*UU(:,1));
    %Del3(i,j)=U0*UU(2,:)*(UU(3,:)'.*FEV(:));

    Density(1,i,j)=UU(1,:)*(UU(1,:)'.*FEV(:));
    %Density(1,i,j)=UU(4,:)*(UU(4,:)'.*(1-FEV(:)));
    Density(2,i,j)=UU(2,:)*(UU(2,:)'.*FEV(:));
    
    Pol(1,i,j)=2*real( UU(2,:)*(UU(1,:)'.*FEV(:))  );
    Pol(2,i,j)=2*imag( UU(2,:)*(UU(1,:)'.*FEV(:))  );
    Pol(3,i,j)=Density(1,i,j)-Density(2,i,j);
    
    %topology(i,j) ; 
    
    end
end


DiffDel=sum(sum(abs(Del-Delold))/(Nx*Ny));





