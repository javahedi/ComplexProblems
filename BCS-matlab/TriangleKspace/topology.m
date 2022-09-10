function topology(i,j)

global t  kx ky UU ED
global alpha GGk


    dTx(1,1)=2*t*sin(kx(i));
    dTx(1,2)=-2*1i*alpha*cos(kx(i));
    dTx(2,1)= 2*1i*alpha*cos(kx(i));
    dTx(2,2)=2*t*sin(kx(i));

    dTy(1,1)=2*t*sin(ky(j));
    dTy(1,2)=2*alpha*cos(ky(j));
    dTy(2,1)=2*alpha*cos(ky(j));
    dTy(2,2)=2*t*sin(ky(j));
    
     dHx=[ dTx zeros(2); zeros(2) dTx.']/2;
    dHy=[ dTy zeros(2); zeros(2) dTy.']/2;
    
    GGk(1,1,i,j)= ( (UU(:,1)'*dHx*UU(:,4))*(UU(:,4)'*dHy*UU(:,1)) ...
                   -(UU(:,1)'*dHy*UU(:,4))*(UU(:,4)'*dHx*UU(:,1)) ) / (ED(1)-ED(4))^2 ...
                 +( (UU(:,1)'*dHx*UU(:,3))*(UU(:,3)'*dHy*UU(:,1)) ...
                   -(UU(:,1)'*dHy*UU(:,3))*(UU(:,3)'*dHx*UU(:,1)) ) / (ED(1)-ED(3))^2;

    GGk(2,1,i,j)= ( (UU(:,2)'*dHx*UU(:,4))*(UU(:,4)'*dHy*UU(:,2)) ...
                   -(UU(:,2)'*dHy*UU(:,4))*(UU(:,4)'*dHx*UU(:,2)) ) / (ED(2)-ED(4))^2 ...
                 +( (UU(:,2)'*dHx*UU(:,3))*(UU(:,3)'*dHy*UU(:,2)) ...
                   -(UU(:,2)'*dHy*UU(:,3))*(UU(:,3)'*dHx*UU(:,2)) ) / (ED(2)-ED(3))^2;

    GGk(1,2,i,j)= ( (UU(:,4)'*dHx*UU(:,1))*(UU(:,1)'*dHy*UU(:,4)) ...
                   -(UU(:,4)'*dHy*UU(:,1))*(UU(:,1)'*dHx*UU(:,4)) ) / (ED(4)-ED(1))^2 ...
                 +( (UU(:,4)'*dHx*UU(:,2))*(UU(:,2)'*dHy*UU(:,4)) ...
                   -(UU(:,4)'*dHy*UU(:,2))*(UU(:,2)'*dHx*UU(:,4)) ) / (ED(4)-ED(2))^2;
            
    GGk(2,2,i,j)= ( (UU(:,3)'*dHx*UU(:,1))*(UU(:,1)'*dHy*UU(:,3)) ...
                   -(UU(:,3)'*dHy*UU(:,1))*(UU(:,1)'*dHx*UU(:,3)) ) / (ED(3)-ED(1))^2 ...
                 +( (UU(:,3)'*dHx*UU(:,2))*(UU(:,2)'*dHy*UU(:,3)) ...
                   -(UU(:,3)'*dHy*UU(:,2))*(UU(:,2)'*dHx*UU(:,3)) ) / (ED(3)-ED(2))^2;

    GGk(:,:,i,j)=1i*GGk(:,:,i,j);
