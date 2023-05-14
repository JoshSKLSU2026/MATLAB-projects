%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Joshua Kruger
% ME 2543--Simulation Methods
% Spring 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1
clear; clc;

%given values
%known constants for water
rho = 62.4; %(pounds/cubic foot)
mu = .01; %(pascal seconds)

D_1 = .1667;
D_2 = .1667; %Diameter Values (ft)
D_3 = .1667;
L_1= 10;
L_2= 10; %Pipe Length Values (ft)
L_3 = 10;
e_1 = .00085;
e_2 = .00085; %Effective surface roughness values (ft)
e_3 = .00085; 

%friction value
f = .043;

%Flow rate (ft^3/s)
Q_big = 1.671 ;


%Cross sectional areas 
A_1 = pi*(D_1/2)^2;
A_2 = pi*(D_2/2)^2;
A_3 = pi*(D_3/2)^2;

%Flow Velocities
% V_1 = Q_1/A_1 ;
% V_2 = Q_2/A_2 ;
% V_3 = Q_3/A_3 ;

%initial value of delta P
deltP= 1 ;

fun2 = @funct;
Q0 = [1,1,1,1];
Q = fsolve(fun2,Q0)
%functions 
Re = functi(Q)
%Friction Values
%Friction1
if Re(1)<=2300
    f_1=64/Re(1)
else
    f_1 = .25*(log(( (e_1/D_1) /3.7) +5.74/ (Re(1)^.9) ))^-2
end

%Friction2
if Re(2)<=2300
    f_2=64/Re(2)
else
    f_2 = .25*(log(( (e_2/D_3) /3.7) +5.74/ (Re(2)^.9) ))^-2
end

%Friction3
if Re(3)<=2300
    f_3=64/Re(3)
else
    f_3 = .25*(log(( (e_3/D_3) /3.7) +5.74/ (Re(3)^.9) ))^-2
end

%PROBLEM ONE EQUATION SOLVER
function g= funct(Q) 

%friction value
syms L_1 L_2 L_3 D_1 D_2 D_3 f rho Q_big A_1 

rho = 62.4; %(pounds/cubic foot)
mu = .01; %(pascal seconds)

D_1 = .1667;
D_2 = .1667; %Diameter Values (ft)
D_3 = .1667;
L_1= 10;
L_2= 10; %Pipe Length Values (ft)
L_3 = 10;
e_1 = .00085;
e_2 = .00085; %Effective surface roughness values (ft)
e_3 = .00085; 

%friction value
f_1 = .043;
f_2 = .043;
f_3 = .043;
%Flow rate (ft^3/s)
Q_big = 1.671 ;

%Cross sectional areas 
A_1 = pi*(D_1/2)^2;
A_2 = pi*(D_2/2)^2;
A_3 = pi*(D_3/2)^2;

g = [f_1*(L_1/D_1)*((Q(1)/A_1)^2/2)-Q(4)/rho;
    f_2*(L_2/D_2)*((Q(2)/A_1)^2/2)-Q(4)/rho; 
    f_3*(L_3/D_3)*((Q(3)/A_1)^2/2)-Q(4)/rho; 
    Q(1)+Q(2)+Q(3)==Q_big;] 
end

%Reynolds Number

function Re = functi(Q)
rho = 62.4; %(pounds/cubic foot)
mu = .01; %(pascal seconds)

D_1 = .1667;
D_2 = .1667; %Diameter Values (ft)
D_3 = .1667;
L_1= 10;
L_2= 10; %Pipe Length Values (ft)
L_3 = 10;
e_1 = .00085;
e_2 = .00085; %Effective surface roughness values (ft)
e_3 = .00085;

%friction value
f = .043;

%Flow rate (ft^3/s)
Q_big = 1.671 ;

%Cross sectional areas
A_1 = pi*(D_1/2)^2;
A_2 = pi*(D_2/2)^2;
A_3 = pi*(D_3/2)^2;



Re(1) = (rho*(Q(1)/A_1)*D_1)/mu;
Re(2) = (rho*(Q(2)/A_2)*D_2)/mu;
Re(3) = (rho*(Q(3)/A_3)*D_3)/mu;
end



%Reynolds number

% Re_1 = (rho*V_1*D_1)/mu ;
% Re_2 = (rho*V_2*D_2)/mu ;
% Re_3 = (rho*V_3*D_3)/mu ;



