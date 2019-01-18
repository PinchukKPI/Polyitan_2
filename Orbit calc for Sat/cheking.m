X  = -4768.6135350850209000;
Y  = -4115.9643385908348000;
Z = -2257.5370510185476;
Vx = -5.24462598;
Vy = 3.74468119;
Vz = 4.24052533;
%| V_u = -5.24462598 V_u = 3.74468119 V_u = 4.24052533  
% x = -4768.6135350850209000 y = -4115.9643385908348000 z = -2257.5370510185476 
 



mu_z = 0.3986*10^15;

mu_z_ = mu_z * 10^(-9.);
%
r = sqrt(X^2. + Y^2. + Z^2.);
V = sqrt(Vx^2. + Vy^2. + Vz^2.);

a = mu_z_*r/(2*mu_z_ - r*V*V);
p = ((Z*Vx - X*Vz)^2. + (X*Vy-Y*Vx)^2. + (Z*Vy-Y*Vz)^2.)/mu_z_;
e = sqrt(1 - p/a);

cos_i = (Z*Vx-X*Vz)/(sqrt(mu_z_*p));
sin_i = sqrt((X*Vy-Y*Vx)^2. + (Z*Vy-Y*Vz)^2.)/(sqrt(mu_z_*p));
i = identification_angel(cos_i, sin_i);

cos_G = (Z*Vy-Y*Vz)/(sin_i*sqrt(mu_z_*p));
sin_G = (X*Vy-Y*Vx)/(sin_i*sqrt(mu_z_*p));
G = identification_angel(cos_G, sin_G);

% cos_u = (Z*cos_G + X*sin_G)/(r*cos_i);
% sin_u = (X*cos_G - Z*sin_G)/(r*cos_i);
cos_u = (Z*cos_G + X*sin_G)/r;
sin_u = Y*(1/sin_i)/r;
u = identification_angel(cos_u, sin_u);

cos_v = (p - r)/(e*r);
sin_v = sqrt(p/mu_z_)*(X*Vx + Y*Vy + Z*Vz)/(e*r);
v = identification_angel(cos_v, sin_v);

wp = u - v;

%orbit_ock = [a, e, i, G, wp, u, v];

fprintf('a = %8.4f e = %8.4f i = %8.4f G = %8.4f wp = %8.4f u = %8.4f v = %8.4f \n', a, e, i, G, wp, u, v);

X  = X * 1000;
Y  = Y * 1000;
Z  = Z * 1000;
Vx = Vx * 1000;
Vy = Vy * 1000;
Vz = Vz * 1000;

mu_z = 0.3986*10^15;
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

%orbit_ock = [a, e, i, G, wp, u, v];

fprintf('a = %8.4f e = %8.4f i = %8.4f G = %8.4f wp = %8.4f u = %8.4f v = %8.4f \n', a, e, i, G, wp, u, v);

