function [vl, vr]=potential(i, it, block)

global Nmid V0

if block==0
  vl=(i+1-Nmid-.5)^2;
  vr=vl;
else
  vl=(i+1-Nmid-.5)^2;
  vr=(i+2-Nmid-.5)^2;
end

if it>0
  vl=(i+1-Nmid-.5)^2;
 if block==0
  vr=vl;
 else
  vr=(i+2-Nmid-.5)^2;
 end
end
  
