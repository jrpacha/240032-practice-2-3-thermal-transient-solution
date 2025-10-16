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
Uini = @(x) 10 + 50*sin(pi*x/2);

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
%Initial Conditions (for the transient solution)
u_ini = Uini(nodes);
%um = u_ini(freeNodes);

%Reamrk: note that in this case in point the essential BC 
%        do not change with t

%Boundary conditions (BC): Natural BC
Q = zeros(numNod, 1); %initialize the global Q vector
Q(freeNodes) = 0; %not necessary, since Q has been initalised to 0.

%Reduced system
Qm = Q(freeNodes) - K(freeNodes,fixedNodes)*uStat(fixedNodes);
u_present = u_ini(freeNodes); %set inital condition

%
%System of equtions: see equation (9)
%
Am = M(freeNodes,freeNodes)+dt*alpha*K(freeNodes,freeNodes);
Bm = M(freeNodes,freeNodes)-dt*(1-alpha)*K(freeNodes,freeNodes);
Fm = dt*Qm;

%Graphical output
    figure()
    plot(nodes, u_ini); %plot the initial state of u_t at t = tini
    title('Temperature Distribution')
    xlabel('x')
    ylabel('temperature')
    hold on
    plot(nodes, uStat, 'g') %Plot stationary solution

%Iterate to find the nodal values of u at each time-tick
for t =tini:dt:tfin-dt %stop one step before, since we compute one step forward
    nextTime = t + dt;
    u_next = Am\(Bm*u_present+Fm); %solve the sistem to find u at the next time-tick
    %plot current curve
        ff = plot(nodes(freeNodes),u_next);
        tt = text(nodes(end)-0.2, 50, ['t = ' num2str(nextTime)]);
        drawnow;
        pause(0.3)
        delete(tt);
        delete(ff);
    %update values
    u_present = u_next;
end
    plot(nodes(freeNodes),u_present);
    text(nodes(end)-0.2, 50, ['t = ' num2str(nextTime)]);
    hold off
%
%Error
%
%Compute the norm of the difference between the stationary and the transient  solution at , i.e., the error :
u_trans = [u_ini(1); u_present; u_ini(end)];

%Precission: compare the difference between the two solutions
errorTemp = norm(uStat-u_trans, inf);
fprintf("error = ||uStat - u_t||_Inf = %e\n", errorTemp)

%Remark. One might increase the final time, tfin, to get an 
%        smaller error.