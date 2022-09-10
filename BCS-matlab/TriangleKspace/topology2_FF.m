function topology2_FF()

global t  kx ky UUk Ek
global alpha GGk Nx Ny QQ

for i=1:Nx
    for j=1:Ny
    dTx(1,1)=2*t*sin(kx(i)+QQ(1));
    dTx(1,2)=-2*1i*alpha*cos(kx(i)+QQ(1));
    dTx(2,1)= 2*1i*alpha*cos(kx(i)+QQ(1));
    dTx(2,2)=2*t*sin(kx(i)+QQ(1));

    dTy(1,1)=2*t*sin(ky(j)+QQ(2));
    dTy(1,2)=2*alpha*cos(ky(j)+QQ(2));
    dTy(2,1)=2*alpha*cos(ky(j)+QQ(2));
    dTy(2,2)=2*t*sin(ky(j)+QQ(2));
    
    dHx=[ dTx zeros(2); zeros(2) dTx.']/2;
    dHy=[ dTy zeros(2); zeros(2) dTy.']/2;
    
    GGk(1,1,i,j)= ( (UUk(:,1,i,j)'*dHx*UUk(:,4,i,j))*(UUk(:,4,i,j)'*dHy*UUk(:,1,i,j)) ...
                   -(UUk(:,1,i,j)'*dHy*UUk(:,4,i,j))*(UUk(:,4,i,j)'*dHx*UUk(:,1,i,j)) ) / (Ek(1,i,j)-Ek(4,i,j))^2 ...
                 +( (UUk(:,1,i,j)'*dHx*UUk(:,3,i,j))*(UUk(:,3,i,j)'*dHy*UUk(:,1,i,j)) ...
                   -(UUk(:,1,i,j)'*dHy*UUk(:,3,i,j))*(UUk(:,3,i,j)'*dHx*UUk(:,1,i,j)) ) / (Ek(1,i,j)-Ek(3,i,j))^2;

    GGk(2,1,i,j)= ( (UUk(:,2,i,j)'*dHx*UUk(:,4,i,j))*(UUk(:,4,i,j)'*dHy*UUk(:,2,i,j)) ...
                   -(UUk(:,2,i,j)'*dHy*UUk(:,4,i,j))*(UUk(:,4,i,j)'*dHx*UUk(:,2,i,j)) ) / (Ek(2,i,j)-Ek(4,i,j))^2 ...
                 +( (UUk(:,2,i,j)'*dHx*UUk(:,3,i,j))*(UUk(:,3,i,j)'*dHy*UUk(:,2,i,j)) ...
                   -(UUk(:,2,i,j)'*dHy*UUk(:,3,i,j))*(UUk(:,3,i,j)'*dHx*UUk(:,2,i,j)) ) / (Ek(2,i,j)-Ek(3,i,j))^2;

    GGk(1,2,i,j)= ( (UUk(:,4,i,j)'*dHx*UUk(:,1,i,j))*(UUk(:,1,i,j)'*dHy*UUk(:,4,i,j)) ...
                   -(UUk(:,4,i,j)'*dHy*UUk(:,1,i,j))*(UUk(:,1,i,j)'*dHx*UUk(:,4,i,j)) ) / (Ek(4,i,j)-Ek(1,i,j))^2 ...
                 +( (UUk(:,4,i,j)'*dHx*UUk(:,2,i,j))*(UUk(:,2,i,j)'*dHy*UUk(:,4,i,j)) ...
                   -(UUk(:,4,i,j)'*dHy*UUk(:,2,i,j))*(UUk(:,2,i,j)'*dHx*UUk(:,4,i,j)) ) / (Ek(4,i,j)-Ek(2,i,j))^2;
            
    GGk(2,2,i,j)= ( (UUk(:,3,i,j)'*dHx*UUk(:,1,i,j))*(UUk(:,1,i,j)'*dHy*UUk(:,3,i,j)) ...
                   -(UUk(:,3,i,j)'*dHy*UUk(:,1,i,j))*(UUk(:,1,i,j)'*dHx*UUk(:,3,i,j)) ) / (Ek(3,i,j)-Ek(1,i,j))^2 ...
                 +( (UUk(:,3,i,j)'*dHx*UUk(:,2,i,j))*(UUk(:,2,i,j)'*dHy*UUk(:,3,i,j)) ...
                   -(UUk(:,3,i,j)'*dHy*UUk(:,2,i,j))*(UUk(:,2,i,j)'*dHx*UUk(:,3,i,j)) ) / (Ek(3,i,j)-Ek(2,i,j))^2;

    GGk(:,:,i,j)=1i*GGk(:,:,i,j);
    end
end
