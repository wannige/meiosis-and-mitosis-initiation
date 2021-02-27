
% This is the main program of the mathematical model of meiosis initiation.
% The current simulation output is meiosis initiation at low nitrogen level
% (N2=0.06).
% other files are needed to be put in the same folder to run this script
% parameters.m:defines the parameters
% eqn1.m : defines the differential equations
clear all;
% simulation time
Tend=30;
% initial conditions of Ime1, Ime2, Cdk1/Cl3 and Rpd3/Sin3
Xo = [0.1 0 0.1 0.6];
% defining some initial values for parameters useful for the while loop
DD=0; % An element of the time matrix ‘tt’
dtt=.1; % An element of the time matrix ‘tt’
j=2; % pointer to be used in the while loop
SOlX=Xo; % SOlX is the final solution matrix
soltime=DD; % soltime is the final time matrix
% solving the ODE at each time step in a while loop, to checking whether
%the protein values are getting negative values. We impose a positive
%basal value for any negative protein value observed.
while(DD<Tend) % starting the while loop
% timespan
tt = DD:dtt:DD+dtt;
% calling the solver
[T X] = ode45(@eqn1,tt,Xo);
% imposing positive basal levels, if negative protein levels observed
if((X(length(T),1))<0),X((length(T)),1)=0.0001;end
if((X(length(T),2))<0),X((length(T)),2)=0.0001;end
if((X(length(T),3))<0),X((length(T)),3)=0.0001;end
if((X(length(T),4))<0),X((length(T)),4)=0.0001;end
% updating the initial condition matrix for the next step
Xo = [X(length(T),1);X(length(T),2);X(length(T),3);X(length(T),4)];
% Adding the solution to the final solution array
SOlX(j,1)=X(length(T),1);
SOlX(j,2)=X(length(T),2);
SOlX(j,3)=X(length(T),3);
SOlX(j,4)=X(length(T),4);
192
DD=DD+dtt; % Updating the element(DD) of the time matrix(tt) for the next
%step
soltime(j)=DD; % Updating the final time matrix
j=j+1; % Updating the pointer for the next step
end % end of the while loop
% Defining the experimental data to be compared with the simulation output
% Ime1 data extracted from experiments
Ime1Prot=[ 1.8 2.9 4.1 5.9 9.15 11.8];
Ime1Pro=[ 1.66/10 3.61/10 6/10 6.18/10 4.6/10 3.8/10 ];
ER=[5.733 7.0428 4.8298 4.025 5.0281 4.5236];
% Ime2 data extracted from experiments
Ime2Prot=[0 5 11.3 13.5 15.6 18.2 20];
Ime2Pro=[0/10 .8/10 3.2/10 4.8/10 4.7/10 4.1/10 3.3/10];
% Rpd3 data extracted from experiments
Rpd3X = [0 1 3 5 12];
Rpd3Y = [9 2.5081 1.1400 1.1400 1.5635]./10;
% Cln3 data extracted from experiments
CLn3IX = [0 0.5054 1.007 1.5036 1.9831 2.990];
CLn3IY = [1 3.1644 3.7505 5.7943 7.1281 10.430];
CLn3DX = [0 0.5139 1.0124 1.5057 1.9919 2.9940];
CLn3DY = [1 1.8188 1.5448 1.4204 1.8195 1.4207];
% plotting the experimental and in-silico protein levels
h1 = figure(1); % Defining a new figure and a header to access the figure
% Setting the position and size of the figure
set(h1, 'Position', [200 200 500 250]);
% plotting the Ime1 experimental data with error bars.
errorbar(Ime1Prot,Ime1Pro,ER/100,'ro');hold on;
plot(Ime2Prot,Ime2Pro,'k*');hold on;
% plotting the scaled Rpd3 data at low nutrients
plot(Rpd3X,(((Rpd3Y./0.9).*0.6)),'bs'); hold on;
% Scaled experimental Cln3 at high nutrients
% plot(CLn3IX,((CLn3IY)./9*(0.14)+.1),'b*'); hold on;
% Scaled experimental Cln3 at low nutrients
plot(CLn3DX,((CLn3DY)./10),'g*'); hold on;
% plotting the in-silico protein levels
plot(soltime,SOlX(:,1),'r','LineWidth',2),hold on;
plot(soltime,SOlX(:,2),'k.','LineWidth',2,'MarkerSize',4),hold on;
plot(soltime,SOlX(:,3),'g-.','LineWidth',2),hold on;
plot(soltime,SOlX(:,4),'b--','LineWidth',2),hold on;
% defining the axis limits for a better representation of the output
axis([0 30 0 0.7])
% defining the axis lables
xlabel('Time (hours)');
ylabel('Relative protein level');