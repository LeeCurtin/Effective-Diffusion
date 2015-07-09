function [ b ] = LoadAssemblerDiscont2( b1,b2 )
%LoadAssemblerDiscont - Creates the load vector for the
%interfaces between all materials including the hydrophobicity of the drug.
b = [b1;b2];

%Couple paste to brain
b(length(b1)+1) = b(length(b1)) + b(length(b1)+1);
b(length(b1)) = 0;

end