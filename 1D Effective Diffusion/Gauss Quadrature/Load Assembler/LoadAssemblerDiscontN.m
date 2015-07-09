function [ b ] = LoadAssemblerDiscontN( b_temp )
%LoadAssemblerDiscont - Creates the load vector for the
%interfaces between all materials including the hydrophobicity of the drug.

b = b_temp{1};

for i = 2:length(b_temp)
    b = [b;b_temp{i}];
end

%Couple all boundaries

for i = 1:length(b_temp)-1
b(length(b_temp{i})+1) = b(length(b_temp{1})) + b(length(b_temp{i})+1);
b(length(b_temp{i})) = 0;
end


end