function warmup(istart,ti,iend)
global  datafile Nop SopL SopR SOL SOR isweep mBt N

tic;
for i=istart:ti-1
   disp(['warm up no truncation : block: ' int2str(i)])
   [opLe,opRe]=enlarge_blocks(i,0);
   [rhoSB,psi,Es]=super_block(i,opLe,opRe,0); % psi is the target state
   for j=1:Nop
     SopL{i+1,j}=opLe{j};
     SopR{i+1,j}=opRe{j};
   end
end

for i=ti:iend
   disp(['warm up: block: ' int2str(i)])
   [opLe,opRe]=enlarge_blocks(i,0);
   [rhoSB,psi,Es]=super_block(i,opLe,opRe,0); % psi is the target state
   Eg(i+1-ti)=Es(1,1);
   OL=den_mat(rhoSB,i,-1);   % forms reduced density matrix and
   update(i,OL,opLe,-1);
   OR=den_mat(rhoSB,i,+1);   % returns the transformation matrix
   update(i,OR,opRe,+1);

  if i==iend
  Spsi{i+1}=update_psi00(i,psi,-1); % wave function transformation
  end
end
 
  data_warmup=[(ti+1:iend+1)' Eg' zeros(size(Eg,2),1)];
  dlmwrite(datafile,data_warmup,'delimiter',' ','precision', 12);
  twrmp=toc;
%  dlmwrite(strcat('data/time_mBt_',int2str(mBt),'_N_',int2str(N),'.dat'), [twrmp isweep], 'delimiter', '\t', '-append');  % warmup time



