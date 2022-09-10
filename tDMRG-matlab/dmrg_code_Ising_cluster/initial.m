 function initial()
 
 global mBt  Jx Jz cluster HH
 global  N ms   Nop
 global datafile isweep tisweep nsweep ntsweep ti  dt
          
 Jx=0;
 Jz=1.;
 cluster=1;
 HH=[0.,0.,0.];
 N=20;           % system size (super-block) must be an even number! 
 Nop=9;    
 mBt=8;          % number of kept states in a block 

 dt=0.1;         % time-step
 nsweep=6;
 ntsweep=0;    % number of tsweep
 isweep=0;
 tisweep=0;

 ctrlsweep=4;    % sweep number to control the convergence  
 convergencE=10^-12;
