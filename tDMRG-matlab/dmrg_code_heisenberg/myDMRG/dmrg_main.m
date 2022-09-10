% main program: Heisenberg spin 1/2 DMRG

tic
clc 
clear all
format long g

global mBt Nmid 
global  N ms  Sops SopL SopR SEs Nop Spsi
global datafile isweep tisweep nsweep ntsweep ti dt

initial();
Nmid=N/2;  % infite-size algorithm will be produced (sub-blocks)
init_mat();

if(mBt<ms) display('mBt<ms'); return; end

ti=floor(log(mBt)/log(ms));     % number of exact sites

if(Nmid>N-ti-3) display('comment out first --> sweep'); return; end

datafile=strcat('data/dmrg_sweep_mBt_',int2str(mBt),'_N_',int2str(N),'.dat');

% WARM UP

warmup(1,ti,Nmid-1); % infinite-size algorithm until superblock size N=2*Nmid

%FINITE_SIZE SWEEPS

for isweep=1:nsweep 
  
  sweep(Nmid  , N-ti-3);   % --->  From middle to rightmost 
  sweep(N-ti-2,   ti+1);   % <---  From rightmost to leftmost
  sweep(ti    , Nmid-1);   % --->  From leftmost to middle
  
end
 
 sweep(Nmid  , N-ti-3); 
 msweep(N-ti-2,    ti);
  
 toc