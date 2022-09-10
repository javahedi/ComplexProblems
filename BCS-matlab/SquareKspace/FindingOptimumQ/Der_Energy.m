function dE=Der_Energy(QQ)

dQQ=0.001;
dE=zeros(1,2);
dE(1,1)=(Kprogram_FF_New(QQ+[dQQ 0.])-Kprogram_FF_New(QQ))/dQQ
dE(1,2)=(Kprogram_FF_New(QQ+[0. dQQ])-Kprogram_FF_New(QQ))/dQQ

%fsolve(@Der_Energy,[0,0])