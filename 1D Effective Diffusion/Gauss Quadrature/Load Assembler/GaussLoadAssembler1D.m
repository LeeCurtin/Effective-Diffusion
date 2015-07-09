function b = GaussLoadAssembler1D(x,f,Neu,Dir,y,w)
%
%GaussLoadAssembler1D - Creates the load vector b using Gaussian quadrature

n = length(x)-1;
b = zeros(n+1,1);

%For the first hat function
Gauss_x = GaussScaling(y,x(1),x(2));
b(1) = sum( f(Gauss_x).*w.*(x(2) - Gauss_x) ); %Normal formulation

if Neu(3) == 1;
b(1) = b(1) - Neu(1); %Added Neumann Condition
end

if Dir(3) == 1;
b(1) = Dir(1); %Dirichlet Condition
end

%For the last hat function
Gauss_xend = GaussScaling(y,x(n),x(n+1));

b(end) = sum( f(Gauss_xend).*w.*(Gauss_xend-x(n)) ); %Normal formulation

if Neu(4) == 1;
b(end) = b(end) + Neu(2); %Added Neumann Condition
end

if Dir(4) == 1;
b(end) = Dir(2); %Dirichlet Condition
end

%For the other hat functions
for i = 2:n
    
    Gauss_x1 = GaussScaling(y,x(i-1),x(i));
    Gauss_x2 = GaussScaling(y,x(i),x(i+1));

    b(i) = sum( f(Gauss_x1).*w.*(Gauss_x1-x(i-1)) +...
                f(Gauss_x2).*w.*(x(i+1)-Gauss_x2) );
end


end