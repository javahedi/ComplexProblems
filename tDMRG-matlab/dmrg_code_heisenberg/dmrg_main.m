% main program: Heisenberg spin 1/2 DMRG
clc 
clear all
format long g

global mBt Nmid 
global  N ms  Sops SopL SopR SEs Nop
global datafile isweep tisweep nsweep ntsweep ti dt

initial();


Nmid=N/2;  % infite-size algorithm will be produced (sub-blocks)


init_mat();

if(mBt<ms) display('mBt<ms'); return; end

ti=floor(log(mBt)/log(ms))     % number of exact sites


if(Nmid>N-ti-3) display('comment out first --> sweep'); return; end

datafile=strcat('data/dmrg_sweep_mBt_',int2str(mBt),'_N_',int2str(N),'.dat');

% WARM UP

warmup(1,ti,Nmid-1); % infinite-size algorithm until superblock size N=2*Nmid

%FINITE_SIZE SWEEPS

for isweep=1:nsweep 
  tic;
  sweep(Nmid  , N-ti-3);   % ---> 
 
  sweep(N-ti-2,   ti+1);   % <---

  sweep(ti    , Nmid-1);   % --->
  
  tswp=toc;
%  dlmwrite(strcat('data/time_mBt_',int2str(mBt),'_N_',int2str(N),'.dat'), [tswp isweep], 'delimiter', '\t', '-append');
 
%  dlmwrite(strcat('data/E_sweep_mBt_',int2str(mBt),'_N_',int2str(N),'.dat'), [SEs{isweep} isweep ],'-append', 'delimiter', '\t', 'precision', 16) 

%  if(isweep>ctrlsweep)         % controls the convergence of the energy
%    DeltaEs=SEs{isweep-ctrlsweep}-SEs{isweep};
%      if(DeltaEs<=convergencE)
%        break
%      end
%  end 
end

return

 sweep(Nmid  , N-ti-3); 
 msweep(N-ti-2,    ti);

 return
 
 % prepare_tsweep();                 % re-calculates operators for the edge of the blocks for the given U, J, V0

  disp('\n TIME EVOLUTION \n')

for tisweep=1:ntsweep
  
  tic;

  tsweep(N-ti-2,   ti+1,    dt/2);   % <---

  tsweep(ti,     N-ti-3,    dt);     % --->

  tsweep(N-ti-2,   ti+1,    dt/2);   % <---

  tsweep(ti,     N-ti-3,      0.);  % --->

  ttswp=toc;
%  dlmwrite(strcat('data/time_mBt_',int2str(mBt),'_N_',int2str(N),'.dat'),[ttswp, tisweep], 'delimiter', '\t', '-append');  % tsweep time

  msweep(N-ti-2,    ti);

end

