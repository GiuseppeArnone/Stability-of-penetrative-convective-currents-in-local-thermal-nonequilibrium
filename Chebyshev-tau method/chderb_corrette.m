function[vector]=chderb_corrette(length,eval)
% Codice preso da Bourne2003
% This subroutine evaluates the derivatives of a Chebyshev polynomial
% at a boundary point of 1 or 0. It returns the vector of
% coefficients that are associated with it in VECTOR.
% VECTOR is the vector returned
% LENGTH is the length of the vector VECTOR
% EVAL is the point at which the derivative is evaluated  -->[0,1]
vector=zeros(1,length);
 if eval==0
     for i=1:4:length-3
     vector(i)=1;
     end
     for i=2:2:length
         vector(i)=0;
     end
     for i=3:4:length-1
         vector(i)=-1;
     end
 elseif eval==1
   for i=1:length
        vector(i)=1;
   end
 elseif eval==-1
    for i=1:length
        vector(i)=eval^(i+1);
    end
 end
end


