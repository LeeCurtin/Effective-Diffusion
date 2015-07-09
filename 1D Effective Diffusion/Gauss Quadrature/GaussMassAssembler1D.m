function [ M ] = GaussMassAssembler1D( x,Dir,y,w )
%StiffnessAssembler1D - Creates the mass matrix A+R using Gaussian
%quadrature

n=length(x)-1;
M=zeros(n+1,n+1);

%Compute A_ii A_i,i+1 and A_i+1,i
for i = 2:n
    h1 = x(i) - x(i-1);
    h2 = x(i+1) - x(i);
    
    Gauss_x1 = GaussScaling(y,x(i-1),x(i));
    Gauss_x2 = GaussScaling(y,x(i),x(i+1));

    M(i,i) = sum( w.*((Gauss_x1-x(i-1)).^2/h1) ) + sum( w.*((x(i+1)-Gauss_x2).^2/h2) );
    M(i,i+1) = sum( w.*((Gauss_x2-x(i)).*(x(i+1)-Gauss_x2))/h2 );
    M(i+1,i) = M(i,i+1);
end

Gauss_x = GaussScaling(y,x(1),x(2));

%Compute A_11
h = x(2) - x(1);
M(1,1) = sum( w.*((x(i+1)-Gauss_x2).^2/h2) ); %Normal formulation

if Dir(3) == 1;
M(1,1) = 0; %For Dirichlet Condition
end

%Compute A_12 and A_21
h = x(2) - x(1);
M(1,2) = sum( w.*((Gauss_x-x(1)).*(x(2)-Gauss_x))/h ); %Normal formulation
M(2,1) = M(1,2);

if Dir(3) == 1;
M(1,2) = 0; %Dirichlet Condition
end


%Compute A_n+1n+1
Gauss_x = GaussScaling(y,x(n),x(n+1));
h = x(n+1) - x(n);

M(n+1,n+1) = sum( w.*((Gauss_x-x(n)).^2/h) ); %Normal formulation

if Dir(4) == 1;
M(n+1,n+1) = 0; %For Dirichlet Condition
M(n+1,n) = 0;
end

end