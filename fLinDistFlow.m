function [pg,qg] = fLinDistFlow(mpc)
%fLinDistFlow: Executes optimal power flow for a radial distribution system


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
   Ar(:,1) = []; %Reduced Incidence matrix,dimension:nline*nline,
    %bus 1 is considered as reference
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
   pinj = Cg*pg - pc; %Real power injection
   qinj = Cg*qg - qc; %Reactive power injection
   Pline = sdpvar(nline,1); %Real line flow
   Qline = sdpvar(nline,1); %Reactive line flow
   v = vref*ones(nbus,1) + 2*R*pinj + 2*X*qinj; %Bus Voltage
   

%Define Constraints
   Constraints=[pinj==A*Pline, qinj==A*Qline , pgmin<=pg(:)<=pgmax , qgmin<=qg(:)<=qgmax , vmin<=v(:)<=vmax , Pline(:) <=2];

%Define Objective
   Objective = c'*pg*mpc.baseMVA;

%Define options
   options=sdpsettings('verbose',1);

%Solution
   sol=optimize(Constraints,Objective,options);
   
%Display results
   Generation=[pg qg]*mpc.baseMVA;Injections=[pinj qinj]*mpc.baseMVA;Flows=[Pline Qline]*mpc.baseMVA;TotalGen=[sum(pg) sum(qg)]*mpc.baseMVA;TotalLoad=[sum(pc) sum(qc)]*mpc.baseMVA;
   disp('  Cost:        ');display(value(Objective));
   disp('  pg:      qg: ');display(value(Generation));
   disp('  Total Generation:');display(value(TotalGen));
   disp('  Total Consumption:');display(value(TotalLoad));
   %disp('  pinj:    qinj:');display(value(Injections));
   disp('  Pline:   Qline:');display(value(Flows));
   disp('  Voltages:');display(value(v));
end

