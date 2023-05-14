%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project #2
% Lexi Garr, Joshua Kruger, Fabian T Chukwuemeka,Tyler Rogers, Nnamdi Dike
% ME 2543--Simulations Methods
% Spring 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

% Problem Number: 1
%initial mass
M0 =.016999;
[theta,sol] = ode45(@func,[pi 3*pi],[1.013e5 0]);


plot(theta,sol(:,1))
xlabel('Theta (rads)')
ylabel('Pressure (Pa)')
title('Pressure vs Theta')

figure()
plot(theta,sol(:,2))
title('Work vs Theta')
xlabel('Theta (rads)')
ylabel('Work (J)')

figure()

volume_g =zeros(length(theta) ,1);
Temp_g =zeros(length(theta),1);
for i = 1:length(theta)
    volume_g(i,:)=V(theta(i,:));
end
plot(theta,volume_g)
xlabel('Theta (rads)')
ylabel('Volume (m^3)')
title('Volume vs Theta')

figure()

% J/kg/k
R = 287;
for i =1:length(theta)
    Temp_g(i,:) =(sol(i,1)*volume_g(i,:))/(M0*R);
end
plot(theta,Temp_g)
xlabel('Theta (rads)')
ylabel('Temperature (K)')
title('Temperature vs Theta')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 2a 
clc;clear;close all;
[theta_2,solu] = ode45(@func_2,[pi 3*pi],[1.013e5 0]);
[theta_3,solut] = ode45(@func_3,[pi 3*pi],[1.013e5 0]);


%Gas constant
% J/kg/k
R = 287;

figure()
% creating zero arrays
volume_g_2 =zeros(length(theta_2),1);
Temp_g_2 =zeros(length(theta_2),1);
volume_g_3 =zeros(length(theta_3),1);
Temp_g_3 =zeros(length(theta_3),1);

%filling zero arrays with values for volume
for i = 1:length(theta_2)
    volume_g_2(i,:)=V_2(theta_2(i,:));
end


for i = 1:length(theta_3)
    volume_g_3(i,:)=V_3(theta_3(i,:));
end

%Filling empty arrays with values for temperature
for i =1:length(theta_2)
    Temp_g_2(i,:) =(solu(i,1)*volume_g_2(i,:))/(mass_2(theta_2(i,:))*R);
end


for i =1:length(theta_3)
    Temp_g_3(i,:) =(solut(i,1)*volume_g_3(i,:))/(mass_3(theta_3(i,:))*R);
end

plot(theta_2,solu(:,1))
hold on
plot(theta_3,solut(:,1))

xlabel('Theta (rads)')
ylabel('Pressure (Pa)')
title('Pressure vs Theta')
legend('Omega=50rad/s','Omega=100rad/s')
figure()

plot(theta_2,solu(:,2))
hold on
plot(theta_3,solut(:,2))
title('Work vs Theta')
xlabel('Theta (rads)')
ylabel('Work (J)')
legend('Omega=50rad/s','Omega=100rad/s')

figure()


plot(theta_2,volume_g_2)
hold on
plot(theta_3,volume_g_3)
xlabel('Theta (rads)')
ylabel('Volume (m^3)')
title('Volume vs Theta')
legend('Omega=50rad/s','Omega=100rad/s')

figure()


plot(theta_2,Temp_g_2)
hold on
plot(theta_3,Temp_g_3)
xlabel('Theta (rads)')
ylabel('Temperature (K)')
title('Temperature vs Theta')
legend('Omega=50rad/s','Omega=100rad/s')

%Problem 1 functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E2 = func(theta,y)
gamma = 1.4;
r = 8.4;
thetas = (3*pi)/2;
thetab = pi;



% meters
l = 012;
S = 008;
b =09;

% cm squared
V0 = 50;

% Kelvin

Tw = 300;

% J/kg
Qin = 2.8e6;

% J/kg/k
R = 287;

%y(1)=Pressure
%y(2)=Work

%mass blow by
he = 0;


%instantaneous surface area
Aw =(4*V(theta))/b;

omega = 1;%rotational velocity
M0 =.016999; %Initial mass

%temperature
T = (y(1)*V(theta))/(M0*R);

%Pressure and work derivatives, respectively
E2 = [-gamma*(y(1)/V(theta))*(V_prime(theta))+(gamma-1)*((M0*Qin)/V(theta))*(x_prime(theta))-((gamma-1)*he*Aw*(T-Tw))/(V(theta)*omega)+((gamma*mass_prime(theta))/M0)*y(1);
    y(1)*V_prime(theta);];


%Volume function
    function Vtheta = V(theta)
        V0 = 50; %cm^3
        sigma = S/(2*l);
        Vtheta = V0*(1 + ( (r-1) / (2*sigma) )*( 1+sigma*(1-cos(theta))-sqrt(1-sigma^2*sin(theta)^2)));
    end

% x derivative function
    function dxdtheta = x_prime(theta)
        if thetas <= theta && theta <= (thetas+thetab)
            %             dxdtheta = -(pi/( 2* thetab))*cos( ( pi*( theta-thetas ) ) / thetab);
            dxdtheta = sin(theta - (3*pi)/2)/2;
        else %if thetas+thetab < theta && theta < 3*pi
            dxdtheta = 0;

        end
    end

%Mass derivative
    function dMdTheta = mass_prime(theta)

        %constants
        omega = 1;%rotational velocity
        M0 = .016999;%initial mass
        C =0 ;%heat transfer to cylinder

        dMdTheta =-( (C*M0) / (exp( (C*theta-pi*C) /omega)*omega) ) ;%( -(C*omega_prime)-(pi*omega_prime) )*M0*exp( (-C/omega)*(theta-pi) );

    end

%Volume derivative
    function dVdTheta = V_prime(theta)
        r = 8.4;
        l = 012;
        S = 008;
        sigma = S/(2*l);
        % dVdTheta = 185*sin(theta) + (185*cos(theta)*sin(theta))/(3*(1 - sin(theta)^2/9)^(1/2));
        dVdTheta = ((2*r*sqrt(1-(sigma^2)*(sin(theta))^2)*sin(theta))+(sigma*r*sin(2*theta))-(2*sqrt(1-(sigma^2)*(sin(theta))^2)*sin(theta)))/(4*sqrt(1-(sigma^2)*(sin(theta))^2));

    end

end





%Volume function
function Vtheta = V(theta)
% meters
l = 012;
S = 008;
r = 8.4;

V0 = 50; %cm^3
sigma = S/(2*l);
Vtheta = V0*(1 + ( (r-1) / (2*sigma) )*( 1+sigma*(1-cos(theta)-sqrt(1-sigma^2*sin(theta)^2))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Problem 2 Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E2_2 = func_2(theta,y)
gamma = 1.4;
r = 8.4;
thetas = (3*pi)/2;
thetab = pi;

% meters
l = 012;
S = 008;
b = 009;

% cm squared
V0 = 50;

% Kelvin

Tw = 300;

% J/kg
Qin = 2.8e6;

% J/kg/k
R = 287;
%y(1)=Pressure
%y(2)=Work

%mass blow by
he = .0050;

%heat transfer to cylinder wall
C = 0.8;

%instantaneous surface area
Aw =(4*V(theta))/b;

%temperature
T = (y(1)*V(theta))/(mass(theta)*R);

omega = 50;%rotational velocity

%Pressure and work derivatives, respectively
E2_2 = [-gamma*(y(1)/V(theta))*(V_prime(theta))+(gamma-1)*((mass(theta)*Qin)/V(theta))*(x_prime(theta))-((gamma-1)*he*Aw*(T-Tw))/(V(theta)*omega)+((gamma*mass_prime(theta))/mass(theta))*y(1);
    y(1)*V_prime(theta);];


%Volume function
    function Vtheta = V(theta)

        V0 = 50; %cm^3
        sigma = S/(2*l);
        Vtheta = V0*(1 + ( (r-1) / (2*sigma) )*( 1+sigma*(1-cos(theta)-sqrt(1-sigma^2*sin(theta)^2))));
    end

% x derivative function
    function dxdtheta = x_prime(theta)
        if thetas <= theta && theta <= (thetas+thetab)
            dxdtheta = sin(theta - (3*pi)/2)/2;
        else %if thetas+thetab < theta && theta < 3*pi
            dxdtheta = 0;

        end
    end



%mass function
    function M = mass(theta)

        M0 = .016999;%initial mass
        C =.8 ;%heat transfer to cylinder
        omega = 50;%rotational velocity
        M = M0*exp((-C/omega)*(theta-pi));
    end

%Mass derivative
    function dMdTheta = mass_prime(theta)

        %constants
        omega = 50;%rotational velocity
        M0 = .016999;%initial mass
        C =.8 ;%heat transfer to cylinder

        dMdTheta =-( (C*M0) / (exp( (C*theta-pi*C) /omega)*omega) ) ;
    end

%Volume derivative
    function dVdTheta = V_prime(theta)
        r = 8.4;
        l = 012;
        S = 008;
        sigma = S/(2*l);
       % dVdTheta = 185*sin(theta) + (185*cos(theta)*sin(theta))/(3*(1 - sin(theta)^2/9)^(1/2));
       dVdTheta = ((2*r*sqrt(1-(sigma^2)*(sin(theta))^2)*sin(theta))+(sigma*r*sin(2*theta))-(2*sqrt(1-(sigma^2)*(sin(theta))^2)*sin(theta)))/(4*sqrt(1-(sigma^2)*(sin(theta))^2));

    end

end





%Volume function
function Vtheta = V_2(theta)
% meters
l = 012;
S = 008;
r = 8.4;
V0 = 50; %cm^3
sigma = S/(2*l);
Vtheta = V0*(1 + ( (r-1) / (2*sigma) )*( 1+sigma*(1-cos(theta)-sqrt(1-sigma^2*sin(theta)^2))));

end

%mass function
function M = mass_2(theta)
C =.8;
omega = 50;
M0 = .016999;%initial mass
M = M0*exp((-C/omega)*(theta-pi));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Problem 2b Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E2_2b = func_3(theta,y)
gamma = 1.4;
r = 8.4;
thetas = (3*pi)/2;
thetab = pi;

% meters
l = 012;
S = 008;
b = 009;

% cm squared
V0 =50;

% Kelvin

Tw = 300;

% J/kg
Qin = 2.8e6;

% J/kg/k
R = 287;
%y(1)=Pressure
%y(2)=Work

%mass blow by
he = .0050;

%heat transfer to cylinder wall
C = 0.8;

%instantaneous surface area
Aw =(4*V(theta))/b;

%temperature
T = (y(1)*V(theta))/(mass(theta)*R);

omega = 100;%rotational velocity

%Pressure and work derivatives, respectively
E2_2b = [-gamma*(y(1)/V(theta))*(V_prime(theta))+(gamma-1)*((mass(theta)*Qin)/V(theta))*(x_prime(theta))-((gamma-1)*he*Aw*(T-Tw))/(V(theta)*omega)+((gamma*mass_prime(theta))/mass(theta))*y(1);
    y(1)*V_prime(theta);];


%Volume function
    function Vtheta = V(theta)

        V0 = 50; %cm^3
        sigma = S/(2*l);
        Vtheta = V0*(1 + ( (r-1) / (2*sigma) )*( 1+sigma*(1-cos(theta)-sqrt(1-sigma^2*sin(theta)^2))));
    end

% x derivative function
    function dxdtheta = x_prime(theta)
        if thetas <= theta && theta <= (thetas+thetab)
            dxdtheta = sin(theta - (3*pi)/2)/2;
        else 
            dxdtheta = 0;

        end
    end


%mass function
    function M = mass(theta)

        M0 = .016999;%initial mass
        C =.8 ;%heat transfer to cylinder
        omega =100;%rotational velocity
        M = M0*exp( (-C/omega) * (theta-pi) );
    end

%Mass derivative
    function dMdTheta = mass_prime(theta)

        %constants
        omega = 100;%rotational velocity

        M0 = .016999;%initial mass
        C =.8 ;%heat transfer to cylinder

        dMdTheta =-( (C*M0) / (exp( (C*theta-pi*C) /omega)*omega) ) ;

    end

%Volume derivative
    function dVdTheta = V_prime(theta)
        r = 8.4;
        l = 012;
        S = 008;
        sigma = S/(2*l);
        %dVdTheta = 185*sin(theta) + (185*cos(theta)*sin(theta))/(3*(1 - sin(theta)^2/9)^(1/2));
        dVdTheta = ((2*r*sqrt(1-(sigma^2)*(sin(theta))^2)*sin(theta))+(sigma*r*sin(2*theta))-(2*sqrt(1-(sigma^2)*(sin(theta))^2)*sin(theta)))/(4*sqrt(1-(sigma^2)*(sin(theta))^2));

    end

end





%Volume function
function Vtheta = V_3(theta)
% meters
l = 012;
S = 008;
r = 8.4;
V0 = 50;
sigma = S/(2*l);
Vtheta = V0*(1 + ( (r-1) / (2*sigma) )*( 1+sigma*(1-cos(theta)-sqrt(1-sigma^2*sin(theta)^2))));
end

%mass function
function M = mass_3(theta)

M0 = .016999;%initial mass
C =.8 ;%heat transfer to cylinder
omega = 100;%rotational velocity
M = M0*exp( (-C/omega) * (theta-pi) );
end