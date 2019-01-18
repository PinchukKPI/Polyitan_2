function SolarCoordinates_ = SolarCoordinates(T);
% Результат в радианах

%Input data
G = T(1);
M = T(2);
d = T(3);
hours = T(4);
min = T(5);
second = T(6);
%
% G = 2012;
% M = 6;
% d = 15;
% hours = 6;
% min = 25;
% second = 42;
% T = [G M d hours min second];
%

% % G = 1978;
% % M = 11;
% % d = 12;
% % hours = 0;
% % min = 0;
% % second = 0;
% % T = [G M d hours min second];
% % Ref = 229+14/60+39/3600 ; % 229 grad 14 min 39 sec

% % G = 1980;
% % M = 7;
% % d = 27;
% % hours = 0;
% % min = 0;
% % second = 0;
% % T = [G M d hours min second];
% % Ref = 124.108828; %grad

% % G = 1978;
% % M = 7;
% % d = 27%26.5;
% % hours = 0;
% % min = 0;
% % second = 0;
% % T = [G M d hours min second];
% % Ref = 123.629912; %grad

JDres = JulianDay(T);
JD = JDres(1);
%
T1900 = [1900 1 0 0 0 0];
JD1900res = JulianDay(T1900);
JD1900 = JD1900res(1);
%
Tu = (JD - JD1900 - 0.5 + (hours*3600 + min*60 + second)/86400)/36525;
%
% Значение аномалии Земли
Mz_ = (99.69098 + 36000.768925*Tu + 0.001*0.38708*Tu*Tu)*pi/180;
%
%
pi_z = ((101 + 13/60 + 15/3600) + (6189.03/3600)*Tu + (1.63/3600)*Tu*Tu)*pi/180;
%
%
e_z = 0.0167504 - 0.0000418*Tu;
%
%
Mz = Mz_ -pi_z;
%
%
n_z = Mz + (2*e_z - 0.25*e_z*e_z*e_z)*sin(Mz) + 1.25*sin(2*Mz)*e_z*e_z;
%
%
n = pi_z + n_z;
%
%
Us = n - pi;
% (Us - 2* pi * fix(Us/(2*pi)))*180/pi - Ref
%
        %fprintf('lambda solar %12.8f \n',Us);
SolarCoordinates_ = Us;
%