function [] = PLOT_D_EFF_2D( u1,p,t,tspan )
%PLOT_CIRCLE_SQUARE_SYSTEM - plots the solution of the TMZ diffusion
%problem

%CREATE SOLUTION MOVIES
% fig = figure(2);
% u = fig.Units;
% fig.Units = 'normalized';
% % fig.Position = [0.3 0.3 0.7 0.7];
% fig.Color = [1 1 1];


F1 = pdeInterpolant(p,t,u1);

xgrid = -3:0.01:3;
ygrid = -3:0.01:3;
[X,Y] = meshgrid(xgrid,ygrid);

uout1 = evaluate(F1,X,Y);
fig = figure(2);
% colormap(cool);
for i = 1:length(tspan)
    u = fig.Units;
    fig.Units = 'normalized';
    fig.Color = [1 1 1];
    colormap(cool);
    Z1 = reshape(uout1(:,i),size(X));
    surf(xgrid,ygrid,Z1);
    axis([-1.5 1.5 -1.5 1.5 -0.01 1.01]);
    title('Diffusion');
    zlabel('Concentration');
    shading interp;
    
    pause(0.001)
%     pause
end

end



