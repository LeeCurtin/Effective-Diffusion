function [ A ] = StiffnessAssemblerDiscontN( A_temp )
%StiffnessAssemblerDiscont - Creates the stiffness matrix for the
%interfaces between all materials

A = A_temp{1};

for j = 2:length(A_temp)
A = blkdiag(A,A_temp{j});
end

%Couple all boundaries

Row_Val = 0; %Declares the row value of interest that will be used to locate where to couple together the stiffness matrices
Col_Val = 0; %Declares the column value equivalent of the row value. 

for i = 1:length(A_temp) - 1
    Row_Val = Row_Val + size(A_temp{i},1);
    Col_Val = Col_Val + size(A_temp{i},2);
    
    for j = 1:Col_Val
        A(Row_Val + 1,j) = A(Row_Val,j);
        A(Row_Val,j) = 0;
    end
    
    A(Row_Val,Col_Val) = -1;
    A(Row_Val,Col_Val + 1) = 1;

end



end



