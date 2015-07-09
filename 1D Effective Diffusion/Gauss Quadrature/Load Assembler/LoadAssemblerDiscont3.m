function [ b ] = LoadAssemblerDiscont3( b1,b2,b3 )
%LoadAssemblerDiscont - Creates the load vector for the
%interfaces between all materials including the hydrophobicity of the drug.
b = [b1;b2;b3];

%Couple paste to water
b(length(b1)+1) = b(length(b1)) + b(length(b1)+1);
b(length(b1)) = 0;

%Couple water to brain
b(length(b1)+length(b2)+1) = b(length(b1)+length(b2)) + b(length(b1)+length(b2)+1);
b(length(b1)+length(b2)) = 0;


% %Couple water to paste
% b2(1) = b1(end) + b2(1);
% b1(end) = 0;
% 
% %Couple paste to water
% b3(1) = b2(end) + b3(1);
% b2(end) = 0;
% 
% %Couple water to brain
% b4(1) = b3(end) + b4(1);
% b3(end) = 0;
% 
% b = [b1;b2;b3;b4];

end