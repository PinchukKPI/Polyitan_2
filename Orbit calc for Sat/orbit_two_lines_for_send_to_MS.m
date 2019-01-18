function orbit_string = orbit_two_lines_for_send_to_MS(L1, L2, in_year, in_month, in_day, in_hours, in_minutes, in_seconds)
%
clc
format long
global r_p r_a i_orb Rz Rz_ekv mu_z w_z w_s eb_ nu_ I I_1 H w_orb J_2 dUT;      % Глобальные переменные для расчета орбиты и параметры спутника
%
Rz = 6371000.;                                      % m         % Средний радиус Земли
Rz_ekv = 6378000.;                                  % m         % Радиус Земли на экваторе
mu_z = 0.3986*10^15;                                % m^3/c^2   % Гравитационная постоянная Земли
% mu_z = 0.3986004418*10^15;                                % m^3/c^2   % Гравитационная постоянная Земли
w_z = 7.292115854937194e-005;                       % c^-1      % Угловая скорость суточного вращения Земли 
w_s = 0.2*10^(-6);                                  % c^-1      % Угловая скорость вращения Земли вокруг Солнца
eb_ = 23.45*pi/180.;                                % rad       % Наклонение плоскости эклиптики к плоскости экватора
nu_ = 11.5*pi/180.;                                 % rad       % Угол наклона диполя МПЗ к оси вращения Земли
J_2 = 1.0826347*10^-3;                              % -         % коэффициент при второй гармонике разложения геопотенциала в ряд Гаусса
%
dUT = 34.*1;
L_mg = 20. * pi/180.;
%
deg_rad = pi()/180; % перевод с град в рад
rad_deg = 180/pi(); % перевод с рад в град
%
test_pytak = 0; % 1 - задан координатами; 2 - задан кеплером
%
% -----------------------------------------------------------------------
% ------------- Input data of orbit -------------------------------------
% -----------------------------------------------------------------------
% 2017   7  10  11  18  34
% L1 = '42732U 98067MM  17191.47122454 +.00013559 +00000-0 +19952-3 0  9998';
% L2 = '42732 051.6397 281.8584 0008157 346.3725 013.7045 15.55737387006960';

% 2017   7  11  10  26   3
% L1 = '42732U 98067MM  17192.43475192 +.00010612 +00000-0 +15761-3 0  9996';
% L2 = '42732 051.6391 277.0438 0008228 349.9734 010.1092 15.55759122007119';

% 2017   7  12   6  28  30 
% L1 = '42732U 98067MM  17193.26979742 +.00009129 +00000-0 +13650-3 0  9994';
% L2 = '42732 051.6388 272.8709 0008281 353.2669 006.8210 15.55777510007242';

% 2017   7  12  14  10  60
% L1 = '42732U 98067MM  17193.59096661 +.00008678 +00000-0 +13008-3 0  9995';
% L2 = '42732 051.6388 271.2660 0008319 354.4537 005.6361 15.55783304007297';

% 2017   7  13  10  13  26
% L1 = '42732U 98067MM  17194.42600090  .00008449  00000-0  12677-3 0  9992';
% L2 = '42732 051.6383 267.0928 0008339 357.6636 002.4315 15.55798101007429';

% 2017   7  13  19  28  25
% L1 = '42732U 98067MM  17194.81139886  .00009408  00000-0  14030-3 0  9996';
% L2 = '42732  51.6386 265.1668 0008383 359.0190   1.0784 15.55806320  7486';

% 2017   7  16  12  13   9
% L1 = '42732U 98067MM  17197.50913694  .00009739  00000-0  14473-3 0  9990';
% L2 = '42732  51.6387 251.6835 0008434   9.0728 351.0414 15.55852950  7904';

% 2017   7  17  12  53   1
% L1 = '42732U 98067MM  17198.53681185 +.00022487 +00000-0 +32387-3 0  9990';
% L2 = '42732 051.6386 246.5467 0008455 012.9031 347.2176 15.55904356008068';

% 1 42732U 98067MM  17200.52785074 +.00009827 +00000-0 +14536-3 0  9994
% 2 42732 051.6384 236.5944 0008611 019.9176 340.2150 15.55967833008374

% 1 42732U 98067MM  17201.49123655 +.00009263 +00000-0 +13734-3 0  9995
% 2 42732 051.6382 231.7786 0008628 023.5396 336.5989 15.55985924008522

% 2017   7  23  21  28  57

%L1 = '1 42732U 98067MM  17204.89510616 +.00011559 +00000-0 +16915-3 0  9994';
%L2 = '2 42732 051.6373 214.7610 0008678 036.0443 324.1132 15.56057851009058';
%L1 = input('Введите первую строку для POLYITAN-2-SAU = ');
%L2 = input('Введите вторую строку для POLYITAN-2-SAU = ');

% -----------------------------------------------------------------------
% ------------- Input data of orbit -------------------------------------
% -----------------------------------------------------------------------
if test_pytak == 1
%     Я взял первую точку из Ваших данных:
%     Date, dmy        09.06.2008
%     Time, hms      07:50:53.646
%     X, km          -5285,009940 
%     Y, km          -4648,166620 
%     Z, km              0,005020    
%     Vx, km/s          -1,028397 
%     Vy, km/s           1,182667  
%     Vz, km/s           7,455440 
% 
%     и проинтегрировал до даты последней точки:
%     Date, dmy        09.06.2008
%     Time, hms     14:22:59.555
% 
%     получил следующий результат:
%     [X Y Z] =  1.0e+006 *[ -3.865135229660671   5.881818835049584   -0.002055075421646 ], m
%     [Vx Vy Vz] =  1.0e+003 *[1.313371552230462   0.855190909416798   7.455571460255504], m/c

% 
%     Ваши данные для последней точки следующие:
%     X, km           -3864,895400
%     Y, km           5882,447370
%     Z, km             0,000000
%     Vx, km/s        1,314666
%     Vy, km/s          0,853377
%     Vz, km/s        7,455183
    
    DateTime = [2008 06 09 07 50 53.646];
    TIME_FOR_ORBIT = [2008 06 09 14 22 59.555];
    XYZ_u_o = [-1.028397 1.182667 7.455440 -5285.009940 -4648.166620 0.005020];
    f_g = [-5285.009940 -4648.166620 0.005020 -1.028397 1.182667 7.455440]*1000;
else
    % --------------------------------------------------------------------------------------
    % Line 1
    L1_array_str = strread(L1, '%s', 'delimiter', ' ');
    % Y_1L = fix(str2double(L1_array_str(3))/1000);
    % D_1L = str2double(L1_array_str(3))-Y_1L*1000;
    Y_1L = fix(str2double(L1_array_str(4))/1000);
    D_1L = str2double(L1_array_str(4))-Y_1L*1000;
    % Line 2
    L2_array_str = strread(L2, '%s', 'delimiter', ' ');
    
    %i_2L = str2double(L2_array_str(2))*deg_rad;        % [deg] - Наклонение орбиты
    %G_2L = str2double(L2_array_str(3))*deg_rad;        % [deg] - Долгота восходящего узла
    %e_2L = str2double(strcat('0.', L2_array_str(4)));  % [-] - Эксцентриситет
    %wp_2L = str2double(L2_array_str(5))*deg_rad;       % [deg]- Аргумент перицентра 
    %M_2L = str2double(L2_array_str(6))*deg_rad;        % [deg] - Средняя аномалия
    %v24h_2L = fix(str2double(L2_array_str(7))*100000000)/100000000;             % [vitok/day] - Частота обращения (оборотов в день) (среднее движение)
    i_2L = str2double(L2_array_str(3))*deg_rad;        % [deg] - Наклонение орбиты
    G_2L = str2double(L2_array_str(4))*deg_rad;        % [deg] - Долгота восходящего узла
    e_2L = str2double(strcat('0.', L2_array_str(5)));  % [-] - Эксцентриситет
    wp_2L = str2double(L2_array_str(6))*deg_rad;       % [deg]- Аргумент перицентра 
    M_2L = str2double(L2_array_str(7))*deg_rad;        % [deg] - Средняя аномалия
    v24h_2L = fix(str2double(L2_array_str(8))*100000000)/100000000;             % [vitok/day] - Частота обращения (оборотов в день) (среднее движение)
    % --------------------------------------------------------------------------------------
    if Y_1L > 56
        G = 1900 + Y_1L;
    else
        G = 2000 + Y_1L;
    end
    t_UT = unixtime([G 1 0 0 0 0])+86400*D_1L;
    Data_UT = unixtime(t_UT);
    M = Data_UT(2);
    d = Data_UT(3);
    hours = Data_UT(4);
    min = Data_UT(5);
    second = Data_UT(6);
    DateTime = [G M d hours min second];
    %
    i_orb = i_2L;
    w0 = 2*pi()/(86400/v24h_2L);    % [rad/sec] - средняя орбитальная угловая скорость
    % w0 = 2*pi()/(86164.090530833/v24h_2L);    % [rad/sec] - средняя орбитальная угловая скорость
    a = (mu_z/(w0*w0))^(1/3);       % [m] - большая полуось
    v = M_2L + (2*e_2L - (e_2L*e_2L*e_2L)/4)*sin(M_2L) + (5*e_2L*e_2L*sin(2*M_2L))/4 + (13*e_2L*e_2L*e_2L*sin(3*M_2L))/12; % [rad] - истинная аномалия
    %
    XYZ_u_o = XYZ_u_(a/1000, e_2L, i_2L, G_2L, wp_2L, wp_2L + v, v);
    %
    % -----------------------------------------------------------------------
    %         TIME_FOR_ORBIT = [2017   7  11  10  24   3; 2017   7  11  10  25   3; 2017   7  11  10  26   3;];
    % TIME_FOR_ORBIT = [2017   7  12   6  26  00; 2017   7  12   6  26  30; 2017   7  12   6  27  00; 2017   7  12   6  28  30;];
    %TIME_FOR_ORBIT = [2017   7  24 19  16  01;];
    TIME_FOR_ORBIT = [ in_year in_month in_day in_hours in_minutes in_seconds];
% -----------------------------------------------------------------------
end
gama_o_calc = sidereal_time(DateTime.', dUT);
L_o = gama_o_calc + L_mg;
gama_o = gama_o_calc;
Ls_o_calc = SolarCoordinates(DateTime);
Ls_o = Ls_o_calc;

%T_ou = T_(2,wp_o + v_o)*T_(3,i_orb)*T_(2,G_o);                  % Вычисление матрицы перехода от ИСК к ОСК
% T_mu = T_(3, nu_)*T_(2, L_o);                                   % Вычисление матрицы перехода от ИСК к МДСК
% T_su_calc = T_(2, Ls_o_calc)*T_(3, eb_);                                  % Вычисление матрицы перехода от ИСК к СЭСК
T_gu = T_(2, gama_o);                                           % Вычисление матрицы перехода от ИСК к ГСК
% T_gu_calc = T_(2, gama_o_calc);                               % Вычисление матрицы перехода от ИСК к ГСК
% -----------------------------------------------------------------------
% if test_pytak == 1
%     XYZ_u = T_gu.'*[f_g(1)*0.001, f_g(2)*0.001, f_g(3)*0.001]';
%     V_u = T_gu.'*[f_g(4)*0.001, f_g(5)*0.001, f_g(6)*0.001]'+w_z*T_dt(2, gama_o).'*[f_g(1)*0.001, f_g(2)*0.001, f_g(3)*0.001]';
%     % Нахожу параметры орбит в ОСК
%     orbit_ock = orbit_ock_(V_u(1), V_u(2), V_u(3), XYZ_u(1), XYZ_u(2), XYZ_u(3));
% else
    xyz_g_o = T_gu*[XYZ_u_o(4),XYZ_u_o(5),XYZ_u_o(6)]';
    orbit_ock = orbit_ock_(XYZ_u_o(1),XYZ_u_o(2),XYZ_u_o(3), XYZ_u_o(4),XYZ_u_o(5),XYZ_u_o(6));
    v_g_o = T_gu*[XYZ_u_o(1),XYZ_u_o(2),XYZ_u_o(3)]'+w_z*T_dt(2, gama_o)*[XYZ_u_o(4),XYZ_u_o(5),XYZ_u_o(6)]';
    f_g = [xyz_g_o' v_g_o'] * 1000;
% end

%fprintf('| Time UTC = %14.2f | Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Xu0 Yuo Zuo [m] = %10.3f %10.3f %10.3f | Vxuo Vyuo Vzuo [m/c] = %10.3f %10.3f %10.3f |\n', unixtime(DateTime), DateTime, [XYZ_u_o(4),XYZ_u_o(5),XYZ_u_o(6)]*1000, [XYZ_u_o(1),XYZ_u_o(2),XYZ_u_o(3)]*1000);
%fprintf('| Time UTC = %14.2f | Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Xg  Yg  Zg  [m] = %10.3f %10.3f %10.3f | Vxg  Vyg  Vzg  [m/c] = %10.3f %10.3f %10.3f |\n', unixtime(DateTime), DateTime, f_g);
%fprintf('| a = %8.4f, e = %9.7f, i = %7.4f, G = %7.4f, wp = %7.4f, u = %8.5f, v = %8.5f\n', orbit_ock(1), orbit_ock(2), orbit_ock(3)*rad_deg, orbit_ock(4)*rad_deg, orbit_ock(5)*rad_deg, orbit_ock(6)*rad_deg, orbit_ock(7)*rad_deg);
 
%
T_b = unixtime(DateTime);
%                      
i_size = size(TIME_FOR_ORBIT);
%
for i_meas = 1:1:i_size(1)
    % 
    T_e = unixtime(TIME_FOR_ORBIT(i_meas,:));
    %
    Rel = 1e-12;
    Y = ode_OSV(@fun_orbit, T_b, T_e, [gama_o f_g].',Rel);
    gama_o = Y(1);
    f_g(1) = Y(2);
    f_g(2) = Y(3);
    f_g(3) = Y(4);
    f_g(4) = Y(5);
    f_g(5) = Y(6);
    f_g(6) = Y(7);
    %
    T_gu = T_(2, gama_o);
    XYZ_u = T_gu.'*[f_g(1)*0.001, f_g(2)*0.001, f_g(3)*0.001]';
    V_u = T_gu.'*[f_g(4)*0.001, f_g(5)*0.001, f_g(6)*0.001]'+w_z*T_dt(2, gama_o).'*[f_g(1)*0.001, f_g(2)*0.001, f_g(3)*0.001]';
    % Нахожу параметры орбит в ОСК

    orbit_ock = orbit_ock_(V_u(1), V_u(2), V_u(3), XYZ_u(1), XYZ_u(2), XYZ_u(3));
    %
	%fprintf('| Time UTC = %14.2f | Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Xu  Yu  Zu  [m] = %10.3f %10.3f %10.3f | Vxu  Vyu  Vzu  [m/c] = %10.3f %10.3f %10.3f |\n', unixtime(TIME_FOR_ORBIT(i_meas,:)), TIME_FOR_ORBIT(i_meas,:), XYZ_u*1000, V_u*1000);
	orbit_string = sprintf('T_orb_0_sec_UTC = %14.2f  \nVx_0  Vy_0  Vz_0  [m/s] = %10.3f  %10.3f  %10.3f\nX_0  Y_0  Z_0  [m] = %10.3f  %10.3f  %10.3f', unixtime(TIME_FOR_ORBIT(i_meas,:)), f_g(4),f_g(5),f_g(6),f_g(1),f_g(2),f_g(3));
    %orbit_string = sprintf('a = %8.4f, e = %9.7f, i = %7.4f, G = %7.4f, wp = %7.4f, u = %8.5f, v = %8.5f\n', orbit_ock(1), orbit_ock(2), orbit_ock(3)*rad_deg, orbit_ock(4)*rad_deg, orbit_ock(5)*rad_deg, orbit_ock(6)*rad_deg, orbit_ock(7)*rad_deg);
    %
    %fprintf(orbit_string);
    T_b = T_e;
    %return afghf;
end    %
PN_settings_T_orb_00_msec_UTC = uint8(0 );%0;
%     uint8_t PN_settings_T_orb_0_msec_UTC;//Время начальной точки орбиты в UTC (результат решения навигационной задачи)
PN_settings_T_orb_0_msec_UTC = uint8(0 );%0;

massive_to_write_int = [PN_settings_T_orb_00_msec_UTC PN_settings_T_orb_0_msec_UTC   ];
massive_to_write_int = typecast(swapbytes(massive_to_write_int) , 'uint16');
  
  %fprintf('%X \n',massive_to_write_int);
%     uint32_t PN_settings_T_orb_0_sec_UTC;//Время начальной точки орбиты в UTC (результат решения навигационной задачи)
PN_settings_T_orb_0_sec_UTC = uint32(unixtime(TIME_FOR_ORBIT(i_meas,:))); %1524124200);
%     float PN_settings_vx_0;              //X Составляющая Вектора скорости (результат решения навигационной задачи)
PN_settings_vx_0 = single(f_g(4)); %5096.97021484375);
%     float PN_settings_vy_0;              //Y Составляющая Вектора скорости (результат решения навигационной задачи)
PN_settings_vy_0 = single(f_g(5)); %3464.14404296875);
%     float PN_settings_vz_0;              //Z Составляющая Вектора скорости (результат решения навигационной задачи)
PN_settings_vz_0 = single(f_g(6)); %-4070.42700195313);
%     float PN_settings_x_0;               //X Координата (результат решения навигационной задачи)
PN_settings_x_0 = single(f_g(1)); %1095664.375);
%     float PN_settings_y_0;               //Y Координата (результат решения навигационной задачи)
PN_settings_y_0 = single(f_g(2)); %4320547.5);
%     float PN_settings_z_0;               //Z Координата (результат решения навигационной задачи)
PN_settings_z_0 = single(f_g(3)); %5054908.5);

PN_settings  = single([ typecast(PN_settings_T_orb_0_sec_UTC , 'single')  PN_settings_vx_0 PN_settings_vy_0 PN_settings_vz_0 PN_settings_x_0 PN_settings_y_0 PN_settings_z_0]);
massive_to_write_PN_settings = massive_to_write_single(PN_settings, 'uint8', 'uint16');
%fprintf('%X ',[massive_to_write_int massive_to_write_PN_settings]);
string_orbit_setting = num2str([massive_to_write_int massive_to_write_PN_settings], '%4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X');
   
orbit_SCF = fopen('orbit_SCF.scf','w');
     
fprintf(orbit_SCF, '<?xml version="1.0" encoding="windows-1251"?>\r\n');
fprintf(orbit_SCF,'<Requests xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="Requests.xsd">\r\n');
fprintf(orbit_SCF,'<MB_req id="ID23" NameRec="1WorbitOpen" Description="открытие записи в настройки подсистемы навигации"><Mb_adr>1</Mb_adr><Mb_func>WRITE_16</Mb_func><Mb_pos>3787</Mb_pos><Mb_num>1</Mb_num><Mb_wr_data>A5A5</Mb_wr_data><Count>1</Count><Timeout_ms>100</Timeout_ms></MB_req><MB_req id="ID25" NameRec="1WorbitPN" Description="запись орбиты в настройки навигации"><Mb_adr>1</Mb_adr><Mb_func>WRITE_16</Mb_func><Mb_pos>3814</Mb_pos><Mb_num>15</Mb_num><Mb_wr_data>');
fprintf(orbit_SCF,string_orbit_setting);
fprintf(orbit_SCF,'</Mb_wr_data><Count>1</Count><Timeout_ms>100</Timeout_ms></MB_req><MB_req id="ID24" NameRec="1WorbitClose" Description="закрытие записи орбиты в настройки навигации"><Mb_adr>1</Mb_adr><Mb_func>WRITE_16</Mb_func><Mb_pos>3787</Mb_pos><Mb_num>1</Mb_num><Mb_wr_data>5A5A</Mb_wr_data><Count>1</Count><Timeout_ms>100</Timeout_ms></MB_req></Requests>\r\n');

fclose(orbit_SCF);

%
end

