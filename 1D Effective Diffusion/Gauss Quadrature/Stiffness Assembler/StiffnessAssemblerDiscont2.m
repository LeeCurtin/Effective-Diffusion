function [ A ] = StiffnessAssemblerDiscont2( A1,A2,h,p )
%StiffnessAssemblerDiscont - Creates the stiffness matrix for the
%interfaces between all materials including the hydrophobicity of the drug.

A = blkdiag(A1,A2);

%Couple paste brain

for i = 1:size(A1,2)
    A(size(A1,1)+1,i) = A(size(A1,1),i);
    A(size(A1,1),i) = 0;
end

A(size(A1,1),size(A1,2)) = -p/h;
A(size(A1,1),size(A1,2)+1) = 1/h;

    
end



