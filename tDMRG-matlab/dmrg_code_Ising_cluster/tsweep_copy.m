function tsweep(istart, iend, dt)
global J U mBt Sops isweep tisweep Spsi N SopL SopR ti Nop ms Id V0 mu 

 if iend>=istart
    block=-1
    [vl,vr]=potential(ti,1,block);
    HBS=expm(-1i*dt*(-J*(kron(SopL{ti,3}, Sops{2})+kron(SopL{ti,2}, Sops{3}))...
                     +(U/2)*(kron(SopL{ti,1}, Sops{6})+kron(SopL{ti,6}, Sops{1}))...
                     -(mu/2)*(kron(SopL{ti,5},Sops{1})+kron(SopL{ti,1},Sops{5}))...
                     +(V0/2)*(kron(SopL{ti,1},vl*Sops{4})+kron(SopL{ti,4},Sops{1}))));
    U1dt=kron(HBS,kron(Sops{1},SopR{N-ti-2,1}));
    Spsi{ti}=U1dt*Spsi{ti};
    bondstart=iend;
  else
    block=+1
%    [vl,vr]=potential(ti-2,1,block);
%    if ti==1
     % HB=expm(-1i*dt*(SopR{ti,Nop}-(U/2)*Sops{6}+(mu/2)*Sops{5}-(1/2)*(V0*vr*Sops{4})));
       HB=expm(-1i*dt*(SopR{ti,Nop}-(U/2)*SopR{ti,6}+(mu/2)*SopR{ti,5}-(V0/2)*SopR{ti,4}));
%    else
     %HB=expm(-1i*dt*(SopR{ti,Nop}-(U/2)*(kron(Sops{5},SopR{ti-1,1}))));
%     HB=expm(-1i*dt*(SopR{ti,Nop}-(U/2)*(SopR{ti,3}*SopR{ti,2})*(SopR{ti,3}*SopR{ti,2}-SopR{ti,1})/2+(mu/2)*kron(Sops{5}, SopR{ti-1,1})-(1/2)*kron(V0*vr*Sops{4},SopR{ti-1,1})));
%    end
    U1dt=kron(kron(SopL{N-ti-2,1}, Sops{1}),kron(Sops{1},HB));
    Spsi{N-ti-2}=U1dt*Spsi{N-ti-2};
    bondstart=istart;
 end

  
for i=istart:-block:iend

    [vl,vr]=potential(i,1,block);

   Hbond=-J*(kron(Sops{3}, Sops{2})+kron(Sops{2}, Sops{3}))...
          +(U/2)*(kron(Sops{6}, Sops{1})+kron(Sops{1}, Sops{6}))...
          -(mu/2)*(kron(Sops{5},Sops{1})+kron(Sops{1},Sops{5}))...
          +(V0/2)*(kron(vl*Sops{4},Sops{1})+kron(Sops{1},vr*Sops{4}));
    Udt=expm(-1i*Hbond*dt);

   
    [opLe,opRe]=enlarge_blocks(i,block);
    %ope=enlarge_blocks(i,block);
    %Ntot=kron(opLe{4},opRe{1})+kron(opLe{1},opRe{4});
    %Ntot=(Ntot+Ntot')/2;

    psi0=Spsi{i};
   
   % tsuper_block(psi0,i,opLe,opRe,block);
    if( mod(i,2) == mod(bondstart,2) )
        disp([int2str(tisweep) '. tsweep, block: ' int2str(i)])
        psi0=kron(SopL{i,1}, kron(Udt, SopR{N-2-i,1}))*psi0;
    end 

    if block==-1
      [O, rho]=den_mat(psi0,i,block); 
    elseif block==+1
      [O, rho]=den_mat(psi0,N-i-2,block);
    end
    %Nexp=psi0'*Ntot*psi0;

    Spsi{i-block}=update_psi00(i, psi0, block);
    %norm(Spsi{i-block})
    if block==-1
      update(i,O,opLe,block);      
    elseif block==+1
      update(N-i-2,O,opRe,block); 
    end  
  % update(-(i-N/2+1)*block+N/2-1,O,ope{1,(block+3)/2},block)

end

if iend>=istart % block=-1;
   [vl,vr]=potential(ti-1,1,block);
   HSB=expm(-1i*dt*(-J*(kron(Sops{3}, SopR{ti,2})+kron(Sops{2}, SopR{ti,3}))...
                  +(U/2)*(kron(Sops{6},SopR{ti,1})+kron(Sops{1}, SopR{ti,6}))...
                  -(mu/2)*(kron(Sops{5},SopR{ti,1})+kron(Sops{1},SopR{ti,5}))...
                  +(V0/2)*(kron(vr*Sops{4},SopR{ti,1})+kron(Sops{1},SopR{ti,4}))));
   U2dt=kron(kron(SopL{N-ti-2,1}, Sops{1}),HSB);
   Spsi{N-ti-2}=U2dt*Spsi{N-ti-2};

 else % block=+1;

    [vl,vr]=potential(ti,1,block);
     Hbond=-J*(kron(Sops{3}, Sops{2})+kron(Sops{2}, Sops{3}))...
          +(U/2)*(kron(Sops{6}, Sops{1})+kron(Sops{1}, Sops{6}))...
          -(mu/2)*(kron(Sops{5},Sops{1})+kron(Sops{1},Sops{5}))...
          +(V0/2)*(kron(vl*Sops{4},Sops{1})+kron(Sops{1},vr*Sops{4}));
    Udt=expm(-1i*Hbond*dt);

    Spsi{ti}=kron(SopL{ti,1},kron(Udt,SopR{N-ti-2,1}))*Spsi{ti};
    psi0=Spsi{ti};
    O=den_mat(psi0,N-ti-2,block);          % wfp is not neccessary block config is the same as the previous step
    update(N-ti-2,O,opRe,block); 
%    if ti==1
     % HB=expm(-1i*dt*(SopL{ti,Nop}-(U/2)*Sops{6}+(mu/2)*Sops{5}-(1/2)*V0*vl*Sops{4}));
       HB=expm(-1i*dt*(SopL{ti,Nop}-(U/2)*SopL{ti,6}+(mu/2)*SopL{ti,5}-(V0/2)*SopL{ti,4}));
%    else
     %HB=expm(-1i*dt*(SopL{ti,Nop}-(U/2)*(kron(SopL{ti-1,1},Sops{5}))));

%      HB=expm(-1i*dt*(SopL{ti,Nop}-(U/2)*(SopL{ti,3}*SopL{ti,2})*(SopL{ti,3}*SopL{ti,2}-SopL{ti,1})/2+(mu/2)*kron(SopL{ti-1,1},Sops{5})-(1/2)*kron(SopL{ti-1,1},V0*vl*Sops{4})));
%    end
    U2dt=kron(kron(HB,Sops{1}),kron(Sops{1},SopR{N-ti-2,1}));
    Spsi{ti}=U2dt*Spsi{ti};
end

% dlmwrite('junk/nexpectation.dat', [tisweep Nexp], 'delimiter', '\t', '-append')

