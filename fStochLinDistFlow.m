function [pg,qg] = fStochLinDistFlow(mpc)
%fStochLinDistFlow: Executes a stochastic optimal power flow for a radial distribution system
%The grid must contain at least two generators.
%Generator at bus no.1 is considered to be the main station.


%Define Constants

%Number of buses & lines & generators
   nbus = size(mpc.bus(:,1));nbus(:,2)=[]; %Bus Number
   nline = size(mpc.branch(:,1));nline(:,2)=[]; %Line Number
   ng  = size(mpc.gen,1); %Total Gen Number
   on = find(mpc.gen(:,8) > 0); %Which generators are on?
   ngen = mpc.gen(on,1); %Which buses are they at?
   ngon = size(on,1); %Active Gen Number
   Cg = sparse(ngen, (1:ngon)', 1, nbus, ng); %Connection matrix element i, j is 1 if bus i contains gen j
   
%A:Incidence Matrix 
   Ainc = makeIncidence(mpc); %Dimension:nline*nbus
   A = Ainc'; %Transposed Incidence matrix,dimension:nbus*nline
   Ar = Ainc;
   Ar(:,1) = []; %Reduced Incidence matrix,dimension:nline*nline,%bus 1 is considered as reference
   F = inv(Ar); %Inverse Incidence matrix
   
%Calculate R and X matrices
   r = diag(mpc.branch(:,3));
   x = diag(mpc.branch(:,4));
   R = F*r*F'; R = [zeros(nline,1) R]; R = [zeros(1,nbus);R];
   X = F*x*F'; X = [zeros(nline,1) X]; X = [zeros(1,nbus);X];
   
%Loads,voltage reference,generator costs,margins
   pc = mpc.bus(:,3)/mpc.baseMVA;
   qc = mpc.bus(:,4)/mpc.baseMVA;
   qgmax = mpc.gen(:,4)/mpc.baseMVA;
   qgmin = mpc.gen(:,5)/mpc.baseMVA;
   pgmax = mpc.gen(:,9)/mpc.baseMVA;
   pgmin = mpc.gen(:,10)/mpc.baseMVA;
   vmax = mpc.bus(:,12);
   vmin = mpc.bus(:,13);
   c = mpc.gencost(:,6);
   vref = 1;
   

%Define Variables
   pg = sdpvar(ngon,1); %Real power generation
   qg = sdpvar(ngon,1); %Reactive power generation
   sg=[pg qg]; %Complex power generation
   sgbound = sdpvar(4,1); %Uncertainty variable
   pinj = Cg*pg - pc; %Real power injection
   qinj = Cg*qg - qc; %Reactive power injection
   Pline = sdpvar(nline,1); %Real line flow
   Qline = sdpvar(nline,1); %Reactive line flow
   v = vref*ones(nbus,1) + 2*R*pinj + 2*X*qinj; %Bus Voltage
   

   
%Define Constraints
   Constraints=[pinj == A*Pline, qinj == A*Qline , pgmin(1) <= sg(1,1) <= pgmax(1), qgmin(1) <= sg(1,2) <= qgmax(1), sgbound(1,1) <= sg(2:ngon,1) <= sgbound(2,1) , sgbound(3,1) <= sg(2:ngon,2) <= sgbound(4,1), vmin <= v(:) <= vmax , Pline <= 2, uncertain(sgbound,'normal',[pgmin(2);pgmax(2);qgmin(2);qgmax(2)],[.001;.01;.01;.01])];
   

%Define Objective
   Objective = c'*pg*mpc.baseMVA;

%Solution
   sol=optimizer(Constraints,Objective,[],sgbound,sg);
   sol_scenarios=sample(sol,1000);
   opt_sol=sol_scenarios([])*mpc.baseMVA;
   
%Display results
   disp('   pg:     qg:  ');display(value(opt_sol));

end

