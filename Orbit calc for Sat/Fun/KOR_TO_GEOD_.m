function KOR_TO_GEOD = KOR_TO_GEOD_(x, y, z)
% Входные параметры
% x, в ГСК по ГОСТ З 51794-2001 [m]
% y, в ГСК по ГОСТ З 51794-2001 [m]
% z, в ГСК по ГОСТ З 51794-2001 [m]
% Выходные параметры
% B - геодезическая широта в [рад]
% L - геодезическая долгота в [рад]
% H - геодезическая высота в [m]

% clc
%x = -1859936.40265821;  y = 5346324.10938850;   z = 4124046.84375836;
% ПЗ-90.02 по ГОСТ Р 51794-2001
alpha_z = 1/298.25784; % сжатие общеземного эллипсоида в [-]
a_z = 6378136; % большая полуось общеземного эллипсоида в [m]
e_z2 = 2*alpha_z - alpha_z*alpha_z;

D = (z*z + x*x)^0.5;
if D == 0
    B = 0.5*pi*y/abs(y);
    L = 0;
    H = y*sin(B) - a_z*(1 - e_z2*(sin(B))^2)^0.5;
else
    La = abs(asin(x/D));
    %fprintf('!!!!KOR_TO_GEOD  La= %12.8f   \n',   La);
    if x >= 0 && z >= 0
        L = La;
    elseif x >= 0 && z <= 0
        L = pi - La;
    elseif x <= 0 && z >= 0
        L = 2*pi - La;
    elseif x <= 0 && z <= 0
        L = pi + La;
    end
    
       % fprintf('!!!!KOR_TO_GEOD  L= %12.8f   \n',   L);
    if y == 0
        B = 0;
        H = D - a_z;
    else
        r = (x*x + y*y + z*z)^0.5;
        c = asin(y/r);
        p = 0.5*e_z2*a_z/r;
        s1 = 0;
        d = 1;
        d_max = 0.000004; %0.0001;
        while d > d_max
            b = c + s1;
            s2 = asin(p*sin(2*b)/((1 - e_z2*(sin(b))^2)^0.5));
            d = abs(s2 - s1);
            s1 = s2;
        end
        B = b;
        H = D*cos(B) + y*sin(B) - a_z*(1 - e_z2*(sin(B))^2)^0.5;
    end
end

% KOR_TO_GEOD_ = [B*180/pi L*180/pi H];
KOR_TO_GEOD = [B L H];

end