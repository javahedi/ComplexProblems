 function initial()
 
 global mBt  J1 J2 D1 D2 HH DM
 global  N ms   Nop
 global datafile isweep tisweep nsweep ntsweep ti  dt
          
 J1=-1.;
 J2=0;
 D1=1.;
 D2=1.;
 DM=0.05;         %  DM interaction
 HH=[0.0,0.,0.];  % hx hy hz
 N=40;           % system size (super-block) must be an even number! 
 Nop=9;    
 mBt=16;          % number of kept states in a block 

 dt=0.1;         % time-step
 nsweep=6;
 ntsweep=0;    % number of tsweep
 isweep=0;
 tisweep=0;

 ctrlsweep=4;    % sweep number to control the convergence  
 convergencE=10^-12;
