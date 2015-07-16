% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[3 2 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-1 1]);
set(ax,'XTick',[ -1.5,...
 -1.3,...
 -1.1000000000000001,...
 -0.89999999999999991,...
 -0.69999999999999996,...
 -0.5,...
 -0.29999999999999982,...
 -0.099999999999999867,...
 0.099999999999999867,...
 0.29999999999999982,...
 0.5,...
 0.69999999999999996,...
 0.89999999999999991,...
 1.1000000000000001,...
 1.3,...
 1.5,...
]);
set(ax,'YTick',[ -1.5,...
 -1.3,...
 -1.1000000000000001,...
 -0.89999999999999991,...
 -0.69999999999999996,...
 -0.5,...
 -0.29999999999999982,...
 -0.099999999999999867,...
 0.099999999999999867,...
 0.29999999999999982,...
 0.5,...
 0.69999999999999996,...
 0.89999999999999991,...
 1.1000000000000001,...
 1.3,...
 1.5,...
]);
pdetool('gridon','on');

% Geometry description:
pderect([-1.3 -1.1000000000000001 0.89999999999999991 0.69999999999999996],'SQ1');
pderect([-0.89999999999999991 -0.69999999999999996 0.89999999999999991 0.69999999999999996],'SQ2');
pderect([-0.5 -0.29999999999999982 0.89999999999999991 0.69999999999999996],'SQ3');
pderect([-0.099999999999999867 0.099999999999999867 0.89999999999999991 0.69999999999999996],'SQ4');
pderect([0.29999999999999982 0.5 0.89999999999999991 0.69999999999999996],'SQ5');
pderect([0.69999999999999996 0.89999999999999991 0.89999999999999991 0.69999999999999996],'SQ6');
pderect([1.1000000000000001 1.3 0.89999999999999991 0.69999999999999996],'SQ7');
pderect([0.89999999999999991 1.1000000000000001 0.5 0.29999999999999982],'SQ8');
pderect([0.5 0.69999999999999996 0.5 0.29999999999999982],'SQ9');
pderect([0.099999999999999867 0.29999999999999982 0.5 0.29999999999999982],'SQ10');
pderect([-0.29999999999999982 -0.099999999999999867 0.5 0.29999999999999982],'SQ11');
pderect([-0.69999999999999996 -0.5 0.5 0.29999999999999982],'SQ12');
pderect([-1.1000000000000001 -0.89999999999999991 0.5 0.29999999999999982],'SQ13');
pderect([-1.3 -1.1000000000000001 0.099999999999999867 -0.099999999999999867],'SQ14');
pderect([-0.89999999999999991 -0.69999999999999996 0.099999999999999867 -0.099999999999999867],'SQ15');
pderect([-0.5 -0.29999999999999982 0.099999999999999867 -0.099999999999999867],'SQ16');
pderect([-0.099999999999999867 0.099999999999999867 0.099999999999999867 -0.099999999999999867],'SQ17');
pderect([0.29999999999999982 0.5 0.099999999999999867 -0.099999999999999867],'SQ18');
pderect([0.69999999999999996 0.89999999999999991 0.099999999999999867 -0.099999999999999867],'SQ19');
pderect([1.1000000000000001 1.3 0.099999999999999867 -0.099999999999999867],'SQ20');
pderect([-1.1000000000000001 -0.89999999999999991 -0.29999999999999982 -0.5],'SQ21');
pderect([-0.69999999999999996 -0.5 -0.29999999999999982 -0.5],'SQ22');
pderect([-0.29999999999999982 -0.099999999999999867 -0.29999999999999982 -0.5],'SQ23');
pderect([0.099999999999999867 0.29999999999999982 -0.29999999999999982 -0.5],'SQ24');
pderect([0.5 0.69999999999999996 -0.29999999999999982 -0.5],'SQ25');
pderect([0.89999999999999991 1.1000000000000001 -0.29999999999999982 -0.5],'SQ26');
pderect([-1.5 1.5 1.0999999999999999 -0.69999999999999996],'R1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','SQ1+SQ2+SQ3+SQ4+SQ5+SQ6+SQ7+SQ8+SQ9+SQ10+SQ11+SQ12+SQ13+SQ14+SQ15+SQ16+SQ17+SQ18+SQ19+SQ20+SQ21+SQ22+SQ23+SQ24+SQ25+SQ26+R1')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');