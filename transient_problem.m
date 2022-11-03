clearvars
close all

format short e

%Solve the evolution problem
% 
% 
% a2 D_t u - D_x ( a_1 D_x u) = 0, u(0) = 10, u(1) = 60
%
%for kc = 2.5, A = 0.1. 
% 
%Compute the solution from t_ini = 0 to t_fin = 0.2 and dt = 0.001. Use 
%500 linear elements and alpha = 1/2 (Crank-Nicolson) value for the
%scheme in eq.4
%
%Take the following initial condition
%
%u_ini(x) = 10 + 50*sin(Pi*x/2)

%Physical constants
kc = 2.5;
A = 0.1;

%Boundary conditions
uA = 10;
uB = 60;

%Initial condition
u_ini = @(x) 10 + 50*sin(pi*x/2);

%FEM: geometry
a = 0; b = 1; numDiv = 500;
h = (b-a)/numDiv;

tini = 0;
tfin = 0.2;
dt = 0.01;
t = tini:dt:tfin; %Time ticks where the (apprximate) solution is computed
alpha = 0.5;      %Crank-Nicolson

a1 = kc*A;
a2 = a1;

nodes = (a:h:b)'; %nodes
elem = [1:numDiv; 2:numDiv+1]'; %connectivity natrix: linear elements

numNod = size(nodes,1);   %find out number of nodes
numElem = size(elem, 1);  %find out number of elements

%Local stiff matrices (linear elements)
Ke = a1*[1, -1; -1, 1]/h; %local stiff matrix (the same for all elements!)
Me = a2*h*[2, 1; 1, 2]/6; %local mass matrix (the same for all elements!)
% Me = a2*h*[1,0;0,1]/2;  %mass lumping
Fe = [0;0];

%Assembly of global matrices
K = zeros(numNod);        %initialise the global stifness matrix 
M = zeros(numNod);        %initialise the global mass matrix
F = zeros(numNod, 1);     %initialise the global vector F

for e = 1:numElem
    rows = [elem(e,1), elem(e,2)];
    cols = rows;
    %Assembly of the stiffness matrices
    K(rows,cols) = K(rows,cols) + Ke;
    M(rows,cols) = M(rows,cols) + Me;
    %Assembly of the F terms (not necessary here)
    %F(rows) = F(rows) + Fe;
end

%Boundary conditions (BC)
fixedNodes = [1, numNod]; %Fixed nodes (nodes where u is fixed)
freeNodes = setdiff(1:numNod,fixedNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Stationary Solution                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Boundary Conditions (BC): Natural BC
Q = zeros(numNod,1);
uStat = zeros(numNod,1);
Q(freeNodes) = 0; %Not necessary, since Q has been initalised to 0.

%Boundary Conditions (BC): Essential BC
uStat(1) = uA; uStat(numNod) = uB;

%Reduced system
Fm = Q(freeNodes) - K(freeNodes,fixedNodes)*uStat(fixedNodes);
Km = K(freeNodes,freeNodes);

%Solve the reduced system
um = Km\Fm;
uStat(freeNodes) = um;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Transient Solution                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = zeros(numNod,1);

%Boundary conditions (BC): Natural BC
Q(freeNodes) = 0; %Not necessary, since Q has been initalised to 0.

%Initial Conditions
u_t = u_ini(nodes);
um = u_t(freeNodes);

%Reamrk: note that that essential BC are the same for each t

%Reduced system
Q_t = Q(freeNodes) - K(freeNodes,fixedNodes)*u_t(fixedNodes);

%System of equtions
Am = M(freeNodes,freeNodes)+dt*alpha*K(freeNodes,freeNodes);
Bm = M(freeNodes,freeNodes)-dt*(1-alpha)*K(freeNodes,freeNodes);
Fm = dt*Q_t;

%Graphical output
figure()
plot(nodes, u_t); %plot the initial state of u_t at t = tini
title('Temperature Distribution')
xlabel('$x$', 'Interpreter','latex')
ylabel('Temperature','Interpreter','latex')
hold on
plot(nodes, uStat, 'g') %Plot stationary soluiton

%Iterate to find the nodal values of u at each time-tick
for t =tini+dt:dt:tfin
    um = Am\(Bm*um+Fm);
    %plot current curve
    ff = plot(nodes(freeNodes),um);
    tt = text(nodes(end)-0.2,50,['t =' num2str(t)]);
    drawnow;
    pause(0.2)
    delete(tt);
    delete(ff);
end
plot(nodes(freeNodes),um);
text(nodes(end)-0.2,50,['t =' num2str(t)]);
hold off
u_t(freeNodes) = um;
%Precission: compare the difference between the two solutions
errorTemp = norm(uStat-u_t, inf);
fprintf("error = ||uStat - u_t||_Inf = %e\n", errorTemp)