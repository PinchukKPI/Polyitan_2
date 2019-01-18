% -----------------------------------------------------------------------
% ---------------- Begin Функция вычисления матрицы поворота ------------
% -----------------------------------------------------------------------
function orbit_ock = orbit_ock_(Vx, Vy, Vz, X, Y, Z) % Begin function orbit_ock_(Vx, Vy, Vz, X, Y, Z) %
% Функция позволяет вычислить параметры орбиты в ОСК используя координаты и скорость в ИСК
% Исходные данные:
% X, Y, Z - [km] координаты в ИСК,
% Vx, Vy, Vz - [км/сек] скорость в ИСК,
% 
mu_z = 0.3986*10^15 * 10^(-9.);                                % m^3/c^2   % Гравитационная постоянная Земли
%
r = sqrt(X^2. + Y^2. + Z^2.);
V = sqrt(Vx^2. + Vy^2. + Vz^2.);

a = mu_z*r/(2*mu_z - r*V*V);
p = ((Z*Vx - X*Vz)^2. + (X*Vy-Y*Vx)^2. + (Z*Vy-Y*Vz)^2.)/mu_z;
e = sqrt(1 - p/a);

cos_i = (Z*Vx-X*Vz)/(sqrt(mu_z*p));
sin_i = sqrt((X*Vy-Y*Vx)^2. + (Z*Vy-Y*Vz)^2.)/(sqrt(mu_z*p));
i = identification_angel(cos_i, sin_i);

cos_G = (Z*Vy-Y*Vz)/(sin_i*sqrt(mu_z*p));
sin_G = (X*Vy-Y*Vx)/(sin_i*sqrt(mu_z*p));
G = identification_angel(cos_G, sin_G);

% cos_u = (Z*cos_G + X*sin_G)/(r*cos_i);
% sin_u = (X*cos_G - Z*sin_G)/(r*cos_i);
cos_u = (Z*cos_G + X*sin_G)/r;
sin_u = Y*(1/sin_i)/r;
u = identification_angel(cos_u, sin_u);

cos_v = (p - r)/(e*r);
sin_v = sqrt(p/mu_z)*(X*Vx + Y*Vy + Z*Vz)/(e*r);
v = identification_angel(cos_v, sin_v);

wp = u - v;
%

orbit_ock = [a, e, i, G, wp, u, v];
%
end                     % End function orbit_ock_(Vx, Vy, Vz, X, Y, Z) %
% -----------------------------------------------------------------------
% ---------------- End Функция вычисления матрицы поворота --------------
% -----------------------------------------------------------------------
