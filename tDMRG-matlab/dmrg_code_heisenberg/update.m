function update(ib,O,ope,block)

global  SopL SopR Nop mBt

% update left block operators to the new basis & truncate Hilbert space & save

if block==-1 

  for j=1:Nop
    SopL{ib+1,j}=O'*ope{j}*O;
  end

elseif block==1

  for j=1:Nop
    SopR{ib+1,j}=O'*ope{j}*O;
  end

else

  stop 'no such block'

end

