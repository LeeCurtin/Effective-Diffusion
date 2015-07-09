function [ A ] = GaussStiffnessAssembler1D_Diff( x,Dir,w,D,y )
%StiffnessAssembler1D - Creates the stiffness matrix A+R using a mid-point
%approximation

n=length(x)-1;
A=zeros(n+1,n+1);

%Compute A_ii A_i,i+1 and A_i+1,i
for i = 2:n
    h1 = x(i) - x(i-1);
    h2 = x(i+1) - x(i);
    
    Gauss_x1 = GaussScaling(y,x(i-1),x(i));
    Gauss_x2 = GaussScaling(y,x(i),x(i+1));
    
    A(i,i) = sum( w.*D(Gauss_x1)*(1/h1) ) + sum( w.*D(Gauss_x2)*(1/h2) );
    A(i,i+1) = sum( w.*D(Gauss_x2)*(-1/h2) );
    A(i+1,i) = A(i,i+1);
end

%Compute A_11
h = x(2) - x(1);
Gauss_x = GaussScaling(y,x(1),x(2));
A(1,1) = sum( w.*D(Gauss_x)*(1/h) ); %Normal formulation

if Dir(3) == 1;
A(1,1) = 1; %For Dirichlet Condition
end

%Compute A_12 and A_21
h = x(2) - x(1);
A(1,2) = sum( w.*D(Gauss_x)*(-1/h) ); %Normal formulation

if Dir(3) == 1;
A(1,2) = 0; %Dirichlet Condition
end

A(2,1) = sum( w.*D(Gauss_x)*(-1/h) );

%Compute A_n+1n+1
h = x(n+1) - x(n);
Gauss_x = GaussScaling(y,x(n),x(n+1));
A(n+1,n+1) = sum( w.*D(Gauss_x)*(1/h) ); %Normal formulation

if Dir(4) == 1;
A(n+1,n+1) = 1; %For Dirichlet Condition
A(n+1,n) = 0;
end

end