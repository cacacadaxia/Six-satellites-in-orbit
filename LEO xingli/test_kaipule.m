



clear all;
clc;

M = 3.5959; % Mean anomaly [rad]
eps = 1E-9; % Tolerance
e = 0.7892; % Eccentricity
En = M;
Ens = En - (En-e*sin(En)- M)/(1 - e*cos(En));
while ( abs(Ens-En) > eps )
    En = Ens;
    Ens = En - (En - e*sin(En) - M)/(1 - e*cos(En));
end
E = Ens; % Eccentric anomaly E [rad].
fprintf('Eccentric anomaly E = %4.4f [rad]\n',E);