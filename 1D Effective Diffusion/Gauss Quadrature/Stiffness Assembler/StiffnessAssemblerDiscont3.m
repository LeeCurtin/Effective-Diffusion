function [ A ] = StiffnessAssemblerDiscont3( A1,A2,A3,h,p )
%StiffnessAssemblerDiscont - Creates the stiffness matrix for the
%interfaces between all materials including the hydrophobicity of the drug.

A = blkdiag(A1,A2,A3);

%Couple paste water

for i = 1:size(A1,2)
    A(size(A1,1)+1,i) = A(size(A1,1),i);
    A(size(A1,1),i) = 0;
end

A(size(A1,1),size(A1,2)) = -1/h;
A(size(A1,1),size(A1,2)+1) = 1/h;


%Couple water to brain
if isempty(A2)       
else
    for i = 1:size(A1,2)+size(A2,2)
        A(size(A1,1)+size(A2,1)+1,i) = A(size(A1,1)+size(A2,1),i);
        A(size(A1,1)+size(A2,1),i) = 0;
    end
    A(size(A1,1)+size(A2,1),size(A1,2)+size(A2,2)) = -p/h;
    A(size(A1,1)+size(A2,1),size(A1,2)+size(A2,2)+1) = 1/h;
end
    
end

