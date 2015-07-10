function [u,X,t,C_0] = D_EFF_1D(h,tau,tend,varargin)
% A function to solve the following reaction diffusion equation:
%C_t = D(x)C_xx 
%h and tau are domain and time step size respectively, with tend being the end time of
%the simulation.
%Here varargin represents n lengths of subdomains to allow for an arbitrary
%number of subdomains

tic %Starts implementation time

%Create Mesh

x0 = 0;
X{1} = x0:h:varargin{1};

X_VAL_START = 0;

for i = 2:length(varargin)
    X_VAL_START = X_VAL_START + varargin{i-1};
    X_VAL_END = X_VAL_START + varargin{i};
    
    X{i} = X_VAL_START:h:X_VAL_END;

end
x = X{1};

for i = 2:length(X)
    x = [x,X{i}];
end
    
% Mesh that is refined around the interfaces between all materials

% x0 = 0;
% L = l1 + l2; %Length of domain
% r = 4; %refinement factor
% x11 = x0              h:      3*l1/4;
% x12 = 3*l1/4+h:       h/r:    l1;
% x21 = l1:             h/r:    l1+l2/4;
% x22 = l1+l2/4+h:      h:      L;

% x1 = [x11,x12];
% x2 = [x21,x22,x23];

% x = [x1,x2];
% end


%Set time interval
t0 = 0;
t = t0:tau:tend;

%Properties
% k = 0.0625; %Rate of splitting of drug and paste
C_0b = 10; %Initial 'bound' drug
C_0f = 0.1; %Initial 'free' drug
C_0 = C_0b + C_0f; %Total drug in system


% %Source term
% g=zeros(1,length(t));
% for i = 1:length(t)
%     g(i) = k*C_0b*exp(-k*t(i))/l1;
% end

%Boundary Conditions
Neu = [0,0,1,1]; %Neumann Boundary Conditions
Dir = [0,0,0,0]; % Dirichlet Boundary Conditions

function [D] = Diffusivity1(x1) 
D1 = 3600*50.e-5;   %Diffusion coefficient (mm^2/h) in wafer for Carmustine
% D1 = 3600*50.e-3;
for i = 1:length(x1)
        D = D1;
end

end

function [D] = Diffusivity2(x2)
D2 = 3600*50.e-5;  %" " " water for Carmustine
% D2 = 3600*50.e-3;
for i = 1:length(x2)
        D = D2;
end

end

function [f] = Sourcef_2(x)
 f = 0;   

% for m = 1:length(x)
%     if (x(m)>=0) && (x(m)<=l1)
%     f = 1;
%     else f = 0;
%     end
% end

end

%Find roots of Legendre polynomial
N=2;
[y,w] = GaussIntegration(N);

%Compute Stiffness Matrix and Load Vector

A_temp = cell(length(varargin),1);
b_temp = cell(length(varargin),1);
M_temp = cell(length(varargin),1);

for i = 1:2:length(varargin);
i;
A_temp{i} = GaussStiffnessAssembler1D_Diff(X{i},Dir,w,@Diffusivity1,y);
A_temp{i+1} = GaussStiffnessAssembler1D_Diff(X{i+1},Dir,w,@Diffusivity2,y);
A_temp

b_temp{1} = GaussLoadAssembler1D(X{1},@Sourcef_2,Neu,Dir,y,w);
b_temp{2} = GaussLoadAssembler1D(X{2},@Sourcef_2,Neu,Dir,y,w);
b_temp{3} = GaussLoadAssembler1D(X{3},@Sourcef_2,Neu,Dir,y,w);
b_temp{4} = GaussLoadAssembler1D(X{4},@Sourcef_2,Neu,Dir,y,w);
b_temp;

M_temp{1} = GaussMassAssembler1D(X{i},Dir,y,w);
M_temp{2} = GaussMassAssembler1D(X{i+1},Dir,y,w);
M_temp{3} = M_temp{1};
M_temp{4} = M_temp{2};
end


% for i = 1:length(varargin)/2;
%     
% A_temp{1} = GaussStiffnessAssembler1D_Diff(X{i},Dir,w,@Diffusivity1,y);
% A_temp{2} = GaussStiffnessAssembler1D_Diff(X{i+1},Dir,w,@Diffusivity2,y);
% A_temp{3} = A_temp{1};
% A_temp{4} = A_temp{2};
% 
% b_temp{1} = GaussLoadAssembler1D(X{i},@Sourcef_2,Neu,Dir,y,w);
% b_temp{2} = GaussLoadAssembler1D(X{i+1},@Sourcef_2,Neu,Dir,y,w);
% b_temp{3} = b_temp{1};
% b_temp{4} = b_temp{2};
% 
% M_temp{1} = GaussMassAssembler1D(X{i},Dir,y,w);
% M_temp{2} = GaussMassAssembler1D(X{i+1},Dir,y,w);
% M_temp{3} = M_temp{1};
% M_temp{4} = M_temp{2};
% end

A = StiffnessAssemblerDiscontN(A_temp);
b = LoadAssemblerDiscontN(b_temp);
M = StiffnessAssemblerDiscontN(M_temp);

%Solve System at each time point

u = zeros(164,length(t));
for i = 1:2:length(varargin)-1
    for j = 1:length(X{i})
        u(j,1) = C_0f;
    end
    for j = 1:length(X{i+1})
        u(j,1) = C_0b;
    end
end
size(M)
size(A)
size(b)
size(u)
theta = 1;
LHS = M/tau + theta*A;

for n = 1:length(t)-1
    RHS = (M/tau-(1-theta)*A)*u(:,n);% + ((1-theta)*g(n)+theta*g(n+1))*b;
    u(:,n+1) = LHS\RHS;
end



toc %Ends implementation time


end


