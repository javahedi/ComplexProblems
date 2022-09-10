 function initial()
 
 global mBt  Jx Jz hx
 global  N ms   Nop
 global datafile isweep tisweep nsweep ntsweep ti  dt
          
 Jz=1.;
 Jx=1;
 hx=0.;
 N=10;           % system size (super-block) must be an even number! 
 Nop=6;    
 mBt=6;          % number of kept states in a block 

 dt=0.1;         % time-step
 nsweep=6;
 ntsweep=0;    % number of tsweep
 isweep=0;
 tisweep=0;

 ctrlsweep=4;    % sweep number to control the convergence  
 convergencE=10^-12;

