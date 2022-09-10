clear all
clc

J0=1;
delta=0.5;

fileID = fopen('system_info.txt','r');
system  = fscanf(fileID,'%f %f %f');
N=system(1);
L=system(2);
alpha=system(3);
fclose(fileID)

site = sort(randperm(L,N))    
%site = [13   30   41   42   50   66   71   73   75   97]
%site = [1:N]

[Sx,Sy,Sz]=spin_spin(N);


H = sparse(0);

for i=1:N
    for j=i+1:N
         J=J0*abs(site(i)-site(j))^(-alpha);
         H = H + J*(cell2mat(Sx(i))*cell2mat(Sx(j))+...
                    cell2mat(Sy(i))*cell2mat(Sy(j)));
    end
end


%{
for i=1:N-1
         H = H + J0*(cell2mat(Sx(i))*cell2mat(Sx(i+1))+...
                     cell2mat(Sy(i))*cell2mat(Sy(i+1))+...
		     delta*cell2mat(Sz(i))*cell2mat(Sz(i+1)));
end
%}

[v,d]=eigs(H,1,"sa");
clear H;
rho = v(:,1)*v(:,1)';
fprintf(' Ground state energy E0/N: %12.8f\n',d(1)/N)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim=2*ones(1,N);            % input for PartialTrace
%sys=[1:1:length(dim)/2];    % input for PartialTrace
% sys  can be very diverse for different purpose
% here we consider half the chain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file1 = fopen('VnEnt.txt','w');
file2 = fopen('Con.txt','w');
file3 = fopen('SySy.txt','w');

cte=1;
for i=1:N
    for j=i+1:N
        d = abs(site(i)-site(j));
    	sysy(cte) = v(:,1)'*cell2mat(Sy(i))*cell2mat(Sy(j))*v(:,1);
    	cte +=1;
    	fprintf(file3,'%i\t%i\t%i\t%12.8f\n',i,j,d,sysy(end));
     end
end

clear Sx; Sy; Sz;

%{
sys=[1]
cte=1;
for i=2:N-1
        tic
	PTR = PartialTrace(rho,sys,dim);
        toc
        %size(PTR);
	VnEnt(cte) = Entropy(PTR);
	sys=[1:1:i]
        cte +=1;
        fprintf(file1,'%12.8f\n',VnEnt(end));
end
%}

%plot(sys,VnEnt)
cte=1;
for i=1:N
       for j=i+1:N
                d = abs(site(i)-site(j));
        	sys=[1:N];
                %i,j
        	sys(i)=[];
        	sys(j-1)=[];
        	PTR = PartialTrace(rho,sys,dim);
        	%size(PTR)
        	Con(cte) = Concurrence(PTR);
        	cte +=1;
        	fprintf(file2,'%i\t%i\t%i\t%12.8f\n',i,j,d,Con(end));
	end
end

fclose(file1);
fclose(file2);

%loglog([1:1:size(Con)],Con)


