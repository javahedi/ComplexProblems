function [Sx,Sy,Sz] = spin_spin(N)

	sx = sparse(0.5* [0 1; 1 0]);
	sy = sparse(0.5*[0 -1i; 1i 0]);
	sz = sparse(0.5*[1 0; 0 -1]);
	id = sparse(eye(2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate operator Sx, Sy, Sz
        Sx = {};
        Sy = {};
        Sz = {};
	for i_site=1:N
    		if i_site==1
       			X=sx;
       			Y=sy;
       			Z=sz;
     		else
       			X=id;
       			Y=id;
       			Z=id;
     		end     
     		for j_site=2:N
        		if j_site==i_site
            			X=kron(X,sx);
       	    			Y=kron(Y,sy);
            			Z=kron(Z,sz);
         		else
  	    			X=kron(X,id);
            			Y=kron(Y,id);
            			Z=kron(Z,id);
         		end
      		end
        Sx{i_site}=X;
        Sy{i_site}=Y;
        Sz{i_site}=Z;
	end
        %  end of generation spin operator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% end of function
