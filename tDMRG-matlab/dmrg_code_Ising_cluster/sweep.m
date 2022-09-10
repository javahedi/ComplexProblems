function sweep(istart,iend) 
global Nmid datafile isweep N SOR SOL ti Spsi Nop SopR SEs U J mu

if iend>=istart
   block=-1
else
   block=+1
end
  
for i=istart:-block:iend
    
    disp([int2str(isweep) '. sweep, block: ' int2str(i)])

    [opLe,opRe]=enlarge_blocks(i,block);

    [rhoSB,psi,Es]=super_block(i,opLe,opRe,block);
    Eg(1+block*(istart-i))=Es(1,1);
    if block ==-1
      [O, rho]=den_mat(rhoSB,i,block);           % find truncation matrix
    elseif block ==+1
      [O, rho]=den_mat(rhoSB,N-i-2,block);
    end
 
    Spsi{i-block}=update_psi00(i,psi,block);
    
   % varN=trace((opLe{5}*rho)^2)-trace(opLe{5}*rho)^2 % variance of particle numbers in each site

    if block==-1
      update(i,O,opLe,block);         % update left block operators
    elseif block==+1
      update(N-i-2,O,opRe,block); 
    end  

    if block==-1 
     if i==Nmid-1
      SEs{isweep}=Es;
     end
    end
end


 data_warmup=[(istart+1:-block:iend+1)' Eg'   isweep*ones(size(Eg,2),1)]; 
 dlmwrite(datafile,data_warmup,'-append','roffset',1,'delimiter',' ','precision', 16);
