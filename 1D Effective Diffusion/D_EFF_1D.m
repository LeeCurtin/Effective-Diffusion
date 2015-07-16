function [u,x,t,C_0] = D_EFF_1D(h,tau,tend,l1,varargin)
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

X{length(varargin)+1} = X_VAL_END:h:X_VAL_END + l1;

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
C_0b = 0.5; %Initial 'bound' drug
C_0f = 0.5; %Initial 'free' drug
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
    
D1 = 3600*50.e-7;   %Diffusion coefficient (mm^2/h) in wafer for Carmustine
D_EFF = 1/((1/(3600*50.e-7) + 1/(3600*50.e-5))/2); %Effective diffusion

for i = 1:length(x1)
%         D = D_EFF;
        D = D1;
end

end

function [D] = Diffusivity2(x2)
    
D2 = 3600*50.e-5;  %" " " water for Carmustine
D_EFF = 1/((1/(3600*50.e-7) + 1/(3600*50.e-5))/2); %Effective diffusion

for i = 1:length(x2)
%         D = D_EFF;
        D = D2;
end

end

function [D] = Diffusivity3(x2)
D2 = 3600*50.e-6;  %" " " water for Carmustine

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

for l = 1:2:length(varargin);

    A_temp{l} = GaussStiffnessAssembler1D_Diff(X{l},Dir,w,@Diffusivity1,y);
    A_temp{l+1} = GaussStiffnessAssembler1D_Diff(X{l+1},Dir,w,@Diffusivity2,y);

    b_temp{l} = GaussLoadAssembler1D(X{l},@Sourcef_2,Neu,Dir,y,w);
    b_temp{l+1} = GaussLoadAssembler1D(X{l+1},@Sourcef_2,Neu,Dir,y,w);

    M_temp{l} = GaussMassAssembler1D(X{l},Dir,y,w);
    M_temp{l+1} = GaussMassAssembler1D(X{l+1},Dir,y,w);

end

A_temp{length(varargin)+1} = GaussStiffnessAssembler1D_Diff(X{end},Dir,w,@Diffusivity3,y);
b_temp{length(varargin)+1} = GaussLoadAssembler1D(X{end},@Sourcef_2,Neu,Dir,y,w);
M_temp{length(varargin)+1} = GaussMassAssembler1D(X{end},Dir,y,w);

A = StiffnessAssemblerDiscontN(A_temp);
b = LoadAssemblerDiscontN(b_temp);
M = StiffnessAssemblerDiscontN(M_temp);

%Solve System at each time point

u = zeros(length(x),length(t));

store_start = [];

for i = 1:2:length(varargin)-1
    
    for j = 1:length(X{i})
        u(j+length(store_start),1) = C_0b;
    end
    
    store_start = [store_start,X{i}];

    for j = 1:length(X{i+1})
        u(j+length(store_start),1) = C_0f;
    end
    
    store_start = [store_start,X{i+1}];

end

% for j = 1:length(X{1})
%     u(j,1) = C_0b;
% end
% for j = length(X{1}):length(X{1}) + length(X{2})
%     u(j,1) = C_0f;
% end
% for j = length(X{1}) + length(X{2}):length(X{1}) + length(X{2}) + length(X{3})
%     u(j,1) = C_0b;
% end
% for j = length(X{1}) + length(X{2}) + length(X{3}):length(X{1}) + length(X{2}) + length(X{3}) + length(X{4})
%     u(j,1) = C_0f;
% end

theta = 1;
LHS = M/tau + theta*A;

for n = 1:length(t)-1
    RHS = (M/tau-(1-theta)*A)*u(:,n);% + ((1-theta)*g(n)+theta*g(n+1))*b;
    u(:,n+1) = LHS\RHS;
end



toc %Ends implementation time


end


