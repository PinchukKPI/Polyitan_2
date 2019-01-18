function orbit_two_lines()
%
clc
format long
global r_p r_a i_orb Rz Rz_ekv mu_z w_z w_s eb_ nu_ I I_1 H w_orb J_2 dUT;      % ���������� ���������� ��� ������� ������ � ��������� ��������
%
FlagB_calc = 1; % 1 - ������������ ������ ������ ���������� ���� �����
T_DNS1 = [1 0 0; 0 1 0; 0 0 1];
T_DNS2 = [1 0 0; 0 -1 0; 0 0 -1]; 
T_DNS3 = [1 0 0; 0 0 -1; 0 1 0];
T_MAG = [-0.01252644 0.99978636 0.01633441; -0.99975786 -0.00309262 -0.02176191; -0.02641699 0.03130639 0.99915776;];
%
Rz = 6371000.;                                      % m         % ������� ������ �����
Rz_ekv = 6378000.;                                  % m         % ������ ����� �� ��������
mu_z = 0.3986*10^15;                                % m^3/c^2   % �������������� ���������� �����
w_z = 7.292115854937194e-005;                       % c^-1      % ������� �������� ��������� �������� ����� 
w_s = 0.2*10^(-6);                                  % c^-1      % ������� �������� �������� ����� ������ ������
eb_ = 23.45*pi/180.;                                % rad       % ���������� ��������� ��������� � ��������� ��������
nu_ = 11.5*pi/180.;                                 % rad       % ���� ������� ������ ��� � ��� �������� �����
J_2 = 1.0826347*10^-3;                              % -         % ����������� ��� ������ ��������� ���������� ������������� � ��� ������
%
dUT = 34.*1;
L_mg = 20. * pi/180.;
%
deg_rad = pi()/180; % ������� � ���� � ���
d_alpha_b_s = 20; % ����, ���� ��������������� ����� ���������� � ��������� ������� ����������� ���������� ������� � ���������� �����������
% 
fid_data1 = fopen('res_MS2.txt','w');
fclose(fid_data1);
%
f_id_inputorbit = fopen('inputorbit.txt','r'); % �������� ����� *.txt
L1 = fgets(f_id_inputorbit);
L2 = fgets(f_id_inputorbit);
i_L1 = 1;
while ischar(L1)
    % Line 1
    L1_array_str = strread(L1, '%s', 'delimiter', ' ');
    Y_1L(i_L1) = fix(str2double(L1_array_str(4))/1000);
    D_1L(i_L1) = str2double(L1_array_str(4))-Y_1L(i_L1)*1000;
    if Y_1L(i_L1) > 56
        G = 1900 + Y_1L(i_L1);
    else
        G = 2000 + Y_1L(i_L1);
    end
    t_UT(i_L1) = unixtime([G 1 0 0 0 0])+86400*D_1L(i_L1);
    
    % Line 2
    L2_array_str = strread(L2, '%s', 'delimiter', ' ');
    i_2L(i_L1) = str2double(L2_array_str(3))*deg_rad;        % [deg] - ���������� ������
    G_2L(i_L1) = str2double(L2_array_str(4))*deg_rad;        % [deg] - ������� ����������� ����
    e_2L(i_L1) = str2double(strcat('0.', L2_array_str(5)));  % [-] - ��������������
    wp_2L(i_L1) = str2double(L2_array_str(6))*deg_rad;       % [deg]- �������� ���������� 
    M_2L(i_L1) = str2double(L2_array_str(7))*deg_rad;        % [deg] - ������� ��������
    v24h_2L(i_L1) = str2double(L2_array_str(8));             % [vitok/day] - ������� ��������� (�������� � ����) (������� ��������)
    
    L1 = fgets(f_id_inputorbit);
    L2 = fgets(f_id_inputorbit);
    
    i_L1 = i_L1 + 1; 
end
fclose(f_id_inputorbit);
% 
f_id_input = fopen('input.txt','r'); % �������� ����� *.txt
names_str = fgets(f_id_input);
names = strread(names_str, '%s', 'delimiter', '	');
i_size = size(names);
for i_names = 1:1:i_size(1)
    if strcmp(names(i_names),'ID1_SUN1_X')
        ID1_SUN1_X = i_names;
    elseif strcmp(names(i_names),'ID1_SUN1_Y')
        ID1_SUN1_Y = i_names;
    elseif strcmp(names(i_names),'ID1_SUN2_X')
        ID1_SUN2_X = i_names;
    elseif strcmp(names(i_names),'ID1_SUN2_Y')
        ID1_SUN2_Y = i_names;
    elseif strcmp(names(i_names),'ID1_SUN3_X')
        ID1_SUN3_X = i_names;
    elseif strcmp(names(i_names),'ID1_SUN3_Y')
        ID1_SUN3_Y = i_names;
    elseif strcmp(names(i_names),'ID1_ADC_RTC_S')
        ID1_ADC_RTC_S = i_names;
    elseif strcmp(names(i_names),'ID1_ADC_RTC_SS')
        ID1_ADC_RTC_SS = i_names;
    elseif strcmp(names(i_names),'ID1_MagTl_X')
        ID1_MagTl_X = i_names;
    elseif strcmp(names(i_names),'ID1_MagTl_Y')
        ID1_MagTl_Y = i_names;
    elseif strcmp(names(i_names),'ID1_MagTl_Z')
        ID1_MagTl_Z = i_names;
    end
end
Meag_L_str = fgets(f_id_input);
i_L1 = 1;
while ischar(Meag_L_str)
    Meag_L = strread(Meag_L_str, '%s', 'delimiter', '	');
    
    DateTime_meg = unixtime(str2double(Meag_L(ID1_ADC_RTC_S)));
    Meag(i_L1,1:6) = [DateTime_meg(1) DateTime_meg(2) DateTime_meg(3) DateTime_meg(4) DateTime_meg(5) DateTime_meg(6)]; %!!!!!!!!!!!!!!!!!!!!!
        
    Meag(i_L1,7) = str2double(Meag_L(ID1_MagTl_X));
    Meag(i_L1,8) = str2double(Meag_L(ID1_MagTl_Y));
    Meag(i_L1,9) = str2double(Meag_L(ID1_MagTl_Z));
    Meag(i_L1,10) = str2double(Meag_L(ID1_SUN1_X));
    Meag(i_L1,11) = str2double(Meag_L(ID1_SUN1_Y));
    Meag(i_L1,12) = str2double(Meag_L(ID1_SUN2_X));
    Meag(i_L1,13) = str2double(Meag_L(ID1_SUN2_Y));
    Meag(i_L1,14) = str2double(Meag_L(ID1_SUN3_X));
    Meag(i_L1,15) = str2double(Meag_L(ID1_SUN3_Y));
    
    Meag_L_str = fgets(f_id_input);
    i_L1 = i_L1 + 1; 
end
fclose(f_id_input);
%
i_size = size(Meag);
%
T_co_LAOO_before = zeros(3,3); %%%%%%%!!!!!!
T_b = 0;
%
for i_meas = 1:1:i_size(1)
    % 
    % -----------------------------------------------------------------------
    % ------------- Input data of Sun and Mag -------------------------------
    % -----------------------------------------------------------------------
    DateTime_meg = [Meag(i_meas, 1) Meag(i_meas, 2) Meag(i_meas, 3) Meag(i_meas, 4) Meag(i_meas, 5) Meag(i_meas, 6)]; %!!!!!!!!!!!!!!!!!!!!!
    T_e = unixtime(DateTime_meg);
    SUN_X1 = Meag(i_meas, 10)*deg_rad;            % [deg]
    SUN_Y1 = Meag(i_meas, 11)*deg_rad;           % [deg]
    SUN_X2 = Meag(i_meas, 12)*deg_rad;            % [deg]
    SUN_Y2 = Meag(i_meas, 13)*deg_rad;           % [deg]
    SUN_X3 = Meag(i_meas, 14)*deg_rad;            % [deg]
    SUN_Y3 = Meag(i_meas, 15)*deg_rad;           % [deg]
    MAG_X = Meag(i_meas, 7);                 % [nT]
    MAG_Y = Meag(i_meas, 8);                  % [nT]
    MAG_Z = Meag(i_meas, 9);                 % [nT]
    % -----------------------------------------------------------------------
    % ------------- Input data of Sun and Mag -------------------------------
    % -----------------------------------------------------------------------
    %
    B_meas = (T_MAG*([MAG_X MAG_Y MAG_Z]*10^(-9))')';
    %
    if (T_e-T_b) > 0 && (T_e-T_b) <= 5400
        rep = 1;
    else
        rep = 0;
        i_orbit = 1;
        while t_UT(i_orbit) > T_e
            i_orbit = i_orbit + 1;
        end   
        
        % -----------------------------------------------------------------------
        % ------------- Input data of orbit -------------------------------------
        % -----------------------------------------------------------------------
    end
   
    if rep == 0
        Data_UT = unixtime(t_UT(i_orbit));
        M = Data_UT(2);
        d = Data_UT(3);
        hours = Data_UT(4);
        min = Data_UT(5);
        second = Data_UT(6);
        DateTime = [G M d hours min second];
%
        i_orb = i_2L(i_orbit);
        w0 = 2*pi()/(86400/v24h_2L(i_orbit));    % [rad/sec] - ������� ����������� ������� ��������
        a = (mu_z/(w0*w0))^(1/3);       % [m] - ������� �������
        v = M_2L(i_orbit) + (2*e_2L(i_orbit) - (e_2L(i_orbit)*e_2L(i_orbit)*e_2L(i_orbit))/4)*sin(M_2L(i_orbit)) + (5*e_2L(i_orbit)*e_2L(i_orbit)*sin(2*M_2L(i_orbit)))/4 + (13*e_2L(i_orbit)*e_2L(i_orbit)*e_2L(i_orbit)*sin(3*M_2L(i_orbit)))/12; % [rad] - �������� ��������
%
        gama_o_calc = sidereal_time(DateTime.', dUT);
        L_o = gama_o_calc + L_mg;
        gama_o = gama_o_calc;
        Ls_o_calc = SolarCoordinates(DateTime);
        Ls_o = Ls_o_calc;

        %T_ou = T_(2,wp_o + v_o)*T_(3,i_orb)*T_(2,G_o);                  % ���������� ������� �������� �� ��� � ���
        % T_mu = T_(3, nu_)*T_(2, L_o);                                   % ���������� ������� �������� �� ��� � ����
        T_su = T_(2, Ls_o)*T_(3, eb_);                                  % ���������� ������� �������� �� ��� � ����
        % T_su_calc = T_(2, Ls_o_calc)*T_(3, eb_);                                  % ���������� ������� �������� �� ��� � ����
        T_gu = T_(2, gama_o);                                           % ���������� ������� �������� �� ��� � ���
        % T_gu_calc = T_(2, gama_o_calc);                               % ���������� ������� �������� �� ��� � ���
        % -----------------------------------------------------------------------
        XYZ_u_o = XYZ_u_(a/1000, e_2L(i_orbit), i_2L(i_orbit), G_2L(i_orbit), wp_2L(i_orbit), wp_2L(i_orbit) + v, v);
        xyz_g_o = T_gu*[XYZ_u_o(4),XYZ_u_o(5),XYZ_u_o(6)]';
        v_g_o = T_gu*[XYZ_u_o(1),XYZ_u_o(2),XYZ_u_o(3)]'+w_z*T_dt(2, gama_o)*[XYZ_u_o(4),XYZ_u_o(5),XYZ_u_o(6)]';
        f_g = [xyz_g_o' v_g_o'] * 1000;
%
        T_b = unixtime(DateTime);
        % T_e = T_b + 100; 
    end

    Rel = 1e-12;
    Y = ode_OSV(@fun_orbit, T_b, T_e, [gama_o f_g].',Rel);
    T_b = T_e;
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
    % ������ ��������� ����� � ���
    orbit_ock = orbit_ock_(V_u(1), V_u(2), V_u(3), XYZ_u(1), XYZ_u(2), XYZ_u(3));
    %
    wp_o = orbit_ock(5); %Y(1);
    G_o = orbit_ock(4); %Y(2);
    v_o = orbit_ock(7); %Y(3);
    % 
    % ----------------------------------------------------------------------------------------------------------
    % ������ ������� ��������� �������� ���������� ���� �����
    % ----------------------------------------------------------------------------------------------------------
    wp_ = wp_o;
    G_ = G_o;
    v_ = v_o;

    u_ = wp_ + v_; 
    T_ou = T_(2,u_)*T_(3,i_orb)*T_(2,G_);                       %           % ������� �������� �� ��� � ���
    r_sat = sqrt(f_g(1)*f_g(1) + f_g(2)*f_g(2) + f_g(3)*f_g(3));   % m         % ������� ������ ������
    v_sat = sqrt(f_g(4)*f_g(4) + f_g(5)*f_g(5) + f_g(6)*f_g(6));   % m/���     % ������� �������� �� �� ������
    
    if FlagB_calc == 1
        KOR_TO_GEOD = KOR_TO_GEOD_(f_g(1), f_g(2), f_g(3));
        B_GEOD = KOR_TO_GEOD(1); % B
        L_GEOD = KOR_TO_GEOD(2); % L
        H_GEOD = KOR_TO_GEOD(3)*0.001; % H
        T_pg = [-cos(L_GEOD) 0 sin(L_GEOD); cos(B_GEOD)*sin(L_GEOD) sin(B_GEOD) cos(B_GEOD)*cos(L_GEOD); -sin(B_GEOD)*sin(L_GEOD) cos(B_GEOD) -sin(B_GEOD)*cos(L_GEOD)];
        B_GOST25645_126_85 = B_GOST25645_126_85_(B_GEOD, L_GEOD, H_GEOD, DateTime, 2010);
        B_GOST_p = [-B_GOST25645_126_85(2) -B_GOST25645_126_85(3) B_GOST25645_126_85(1)]'*10^(-9);
        T_res = T_ou * T_gu.' * T_pg.';
        B_o_print = T_ou * T_gu.' * T_pg.'*B_GOST_p;
    else
        T_mu = T_(3, nu_)*T_(2, L_o);                                %           % ������� �������� �� ��� � ����
        T_um = T_mu.'; 
        B_o_print = B_o(r_, T_ou, T_um);
    end
    % ----------------------------------------------------------------------------------------------------------
    % ������ ������� ��������� �������� ���������� ���� �����
    % ----------------------------------------------------------------------------------------------------------
    T_uo = T_ou.';
%     q0_o = q(1);
%     q1_o = q(2);
%     q2_o = q(3);
%     q3_o = q(4);
%     T_co = T_qvant(q0_o, q1_o, q2_o, q3_o);                                             % ������� Lamda_T(T_)

%     Ls_ = Ls_o_calc;
%     T_su = T_(2, Ls_)*T_(3, eb_);                               %           % ������� �������� �� ��� � ����
%     
    T_so = T_su*T_uo;    
%     T_cs = T_co*(T_so.');                                                   %           % ������� �������� �� ���� � ���
%     es_c = T_cs(:,3);
%     S_c = es_c;
    T_os = T_so.';
    S_o_print = T_os(:,3);
%     B_c = T_co*B_o_print;
%   
    % -----------------------------------------------
    B_o_calc = B_o_print;
    S_o_calc = S_o_print;
    % -----------------------------------------------
    IsMeas_MDA = norm(B_meas);
    IsCalc_S_calc = 1;
    IsCalc_B_calc = norm(B_o_calc);
    %
    S_meas = [];

    grad_S_1_MDA = NaN;
    IsMeas_S_1_MDA = 0;
    if abs(SUN_X1) < 3*pi() && abs(SUN_Y1) < 3*pi()
        S_1 = T_DNS1*(S_alfa_beta(SUN_X1, SUN_Y1))';
        
        IsMeas_S_1 = norm(S_1);
        grad_S_o_calc_B_o_calc = acosd(S_o_calc'*B_o_calc/(IsCalc_S_calc*IsCalc_B_calc));
        grad_S_1_MDA = acosd(S_1'*B_meas'/(IsMeas_S_1*IsMeas_MDA));
        if abs(grad_S_1_MDA - grad_S_o_calc_B_o_calc) <= d_alpha_b_s
            IsMeas_S_1_MDA = 1;
        end
    	fprintf('| Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Grad_calc Grad_meas1 dGrad = %7.3f %7.3f %7.3f |\n', DateTime_meg, grad_S_o_calc_B_o_calc, grad_S_1_MDA, (grad_S_o_calc_B_o_calc - grad_S_1_MDA));
    end
    
    grad_S_2_MDA = NaN;
    IsMeas_S_2_MDA = 0;
    if abs(SUN_X2) < 3*pi() && abs(SUN_Y2) < 3*pi()
        S_2 = T_DNS2*(S_alfa_beta(SUN_X2, SUN_Y2))';

        IsMeas_S_2 = norm(S_2);
        grad_S_o_calc_B_o_calc = acosd(S_o_calc'*B_o_calc/(IsCalc_S_calc*IsCalc_B_calc));
        grad_S_2_MDA = acosd(S_2'*B_meas'/(IsMeas_S_2*IsMeas_MDA));
        if abs(grad_S_2_MDA - grad_S_o_calc_B_o_calc) <= d_alpha_b_s
            IsMeas_S_2_MDA = 1;
        end
    	fprintf('| Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Grad_calc Grad_meas2 dGrad = %7.3f %7.3f %7.3f |\n', DateTime_meg, grad_S_o_calc_B_o_calc, grad_S_2_MDA, (grad_S_o_calc_B_o_calc - grad_S_2_MDA));
    end
    
    grad_S_3_MDA = NaN;
    IsMeas_S_3_MDA = 0;
    if abs(SUN_X3) < 3*pi() && abs(SUN_Y3) < 3*pi()
        S_3 = T_DNS3*(S_alfa_beta(SUN_X3, SUN_Y3))';

        IsMeas_S_3 = norm(S_3);
        grad_S_o_calc_B_o_calc = acosd(S_o_calc'*B_o_calc/(IsCalc_S_calc*IsCalc_B_calc));
        grad_S_3_MDA = acosd(S_3'*B_meas'/(IsMeas_S_3*IsMeas_MDA));
        if abs(grad_S_3_MDA - grad_S_o_calc_B_o_calc) <= d_alpha_b_s
            IsMeas_S_3_MDA = 1;
        end
    	fprintf('| Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Grad_calc Grad_meas3 dGrad = %7.3f %7.3f %7.3f |\n', DateTime_meg, grad_S_o_calc_B_o_calc, grad_S_3_MDA, (grad_S_o_calc_B_o_calc - grad_S_3_MDA));
    end
    
	if IsMeas_S_1_MDA == 1
        S_meas = [S_meas S_1'];   
    end 
	if IsMeas_S_2_MDA == 1
        S_meas = [S_meas S_2'];   
    end 
	if IsMeas_S_3_MDA == 1
        S_meas = [S_meas S_3'];   
    end 
%
%     IsMeas_S_2 = norm(S_2);
%     IsMeas_S_3 = norm(S_3);
   
    %-----------------------------------------------------------
    % ��������� ����� ����������� ����������
    %-----------------------------------------------------------
	i_s = size(S_meas);
    i_s_end = i_s(2)/3;
    i_b = size(B_meas);
    i_b_end = i_b(2)/3;
    T_co_LAOO_i = zeros(3,3);
    w_co_c_LAOO_i = zeros(3,3);
    
    if i_s_end >= 1 && i_b_end >= 1
        
	for i_b = 1:1:i_b_end
        for i_s = 1:1:i_s_end
                
            B_c_i = B_meas((i_b*3-2):i_b*3)';
            S_c_i = S_meas((i_s*3-2):i_s*3)';
            T_co_LAOO = LAOO(B_c_i, B_o_print, S_c_i, S_o_print);
            T_co_LAOO_i = T_co_LAOO_i + T_co_LAOO;
            if norm(T_co_LAOO_before) == 0
                T_co_LAOO_before = T_co_LAOO;
                T_control = 1;
            else
                t_i = unixtime([Meag(i_meas, 1) Meag(i_meas, 2) Meag(i_meas, 3) Meag(i_meas, 4) Meag(i_meas, 5) Meag(i_meas, 6)]);
                t_i_1 = unixtime([Meag(i_meas-1, 1) Meag(i_meas-1, 2) Meag(i_meas-1, 3) Meag(i_meas-1, 4) Meag(i_meas-1, 5) Meag(i_meas-1, 6)]);
                T_control = t_i - t_i_1;
            end
            T_co_LAOO_st = (T_co_LAOO - T_co_LAOO_before)/(T_control);
            
            T_co_LAOO_1 = inv(T_co_LAOO);
            w_co_c_LAOO = T_co_LAOO_st*T_co_LAOO_1;
            w_co_c_LAOO_i = w_co_c_LAOO_i + w_co_c_LAOO;
                
        end
    end
        
    if i_s_end ~= 0 && i_b_end ~= 0
        T_co_LAOO = T_co_LAOO_i/(i_b*i_s);
        w_co_c_LAOO = w_co_c_LAOO_i/(i_b*i_s);
    
        w_co_c_LAOO1 = -[-w_co_c_LAOO(2,3) w_co_c_LAOO(1,3) -w_co_c_LAOO(1,2)];
        w_co_c_LAOO2 = -[w_co_c_LAOO(3,2) -w_co_c_LAOO(3,1) w_co_c_LAOO(2,1)];
        w_co_c_LAOO_rez = (w_co_c_LAOO1 + w_co_c_LAOO2)*0.5;
%        p_ = a_*(1-e_^2.);                                                       % m          %��������� �������� ������
%        wp_st = 5.*((Rz_ekv*0.001/a_)^3.5)*(5.*(cos(i_orb))^2.-1.)*(pi/180.)/(24.*3600.); % rad/sec %����������� ��������� �������
%        v_st = sqrt(mu_z*p_*1000)/(r_^2.);
%        wp_st_o = [0. wp_st 0.].';
%        v_st_o = [0. v_st 0.].';
        G_st = -10.*((Rz_ekv*0.001/a)^3.5)*(cos(i_orb))*(pi/180.)/(24.*3600.);           % rad/sec %����������� ����������� ����
        G_st_u = [0. G_st 0.].';

%        u_st_o = wp_st_o + v_st_o;
%        u_st_o = [0. sqrt(mu_z/r_^3.) 0.].';
        u_st_o = [0. (v_sat/r_sat) 0.].';
        %
        w_orb_p = u_st_o; %T_co*[0, norm(V_u)/norm(XYZ_u), 0]';
        %
        w_ou_o = u_st_o + T_ou*G_st_u;                                               % ������ ������� �������� ��� ������������ ���, ���������� � ��������� �� ��� ���
        T_co_LAOO_w_ou_o = T_co_LAOO*w_ou_o;
        w_cu_c_LAOO_rez = w_co_c_LAOO_rez.' + T_co_LAOO_w_ou_o;                             % ������ ������� �������� �������� ������������ ��� � ��������� �� ��� ���

    else
        T_co_LAOO = zeros(3,3);
        w_orb_p = [0 0 0].';
    end
    
    T_co_LAOO_before = T_co_LAOO;
	fprintf('| Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | TETA FI PSI = %7.3f %7.3f %7.3f | w_cu_c_LAOO_rez = %7.3f %7.3f %7.3f |\n', DateTime_meg, ELIR_angle(Lamda_T(T_co_LAOO)), w_cu_c_LAOO_rez*180/pi);

    fid_data1 = fopen('res_MS2.txt','at+');
    fprintf(fid_data1,'| Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | TETA FI PSI [deg] = %7.3f %7.3f %7.3f | w_cu_c_LAOO_rez [deg/sec] = %7.3f %7.3f %7.3f | Grad_calc [deg] = %7.3f | Grad_meas1 dGrad1 [deg] = %7.3f %7.3f | Grad_meas2 dGrad2 [deg] = %7.3f %7.3f | Grad_meas3 dGrad3 [deg] = %7.3f %7.3f |\n', DateTime_meg, ELIR_angle(Lamda_T(T_co_LAOO)), w_cu_c_LAOO_rez*180/pi, grad_S_o_calc_B_o_calc, grad_S_1_MDA, (grad_S_o_calc_B_o_calc - grad_S_1_MDA), grad_S_2_MDA, (grad_S_o_calc_B_o_calc - grad_S_2_MDA), grad_S_3_MDA, (grad_S_o_calc_B_o_calc - grad_S_3_MDA));
    fclose(fid_data1);
    
    else
        fid_data1 = fopen('res_MS2.txt','at+');
        fprintf(fid_data1,'| Y M D h m s = %5.0f %3.0f %3.0f %3.0f %3.0f %3.0f | Measured data is not valid | Grad_calc [deg] = %7.3f | Grad_meas1 dGrad1 [deg] = %7.3f %7.3f | Grad_meas2 dGrad2 [deg] = %7.3f %7.3f | Grad_meas3 dGrad3 [deg] = %7.3f %7.3f |\n', DateTime_meg, grad_S_o_calc_B_o_calc, grad_S_1_MDA, (grad_S_o_calc_B_o_calc - grad_S_1_MDA), grad_S_2_MDA, (grad_S_o_calc_B_o_calc - grad_S_2_MDA), grad_S_3_MDA, (grad_S_o_calc_B_o_calc - grad_S_3_MDA));
        fclose(fid_data1);
    end
        
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    %-----------------------------------------------------------

    end

end


% -----------------------------------------------------------------------
% ---------------- Begin ������� �������� ������������ � ���� ������ ----
% -----------------------------------------------------------------------
function ELIR_angle_ = ELIR_angle(q)                % Begin function ELIR_angle(q) %
% ������� ��������� ����������� � ���� ������, � ������ �������� �������
% ��������� ������������������ ���������: ������, ����, ��������.
% �������� ������: q - ������ �������� ������������.
% �������� �������� - ���� �������, ����� � �������� (ELIR_angle), ����.
del = 180./pi;                                      %       % ��������� �������� � ������ � �������
q0 = q(1);                                          %
q1 = q(2);                                          %
q2 = q(3);                                          %
q3 = q(4);                                          %
T_x = T_qvant(q0, q1, q2, q3);                      %       % ������ ������� �������� � �������������� ������������
TETA_x = del*atan2(T_x(3,1),T_x(3,3));              % grad  % ������ ���� �������
FI_x = del*atan2(-T_x(3,2),sqrt(1.-T_x(3,2)^2.));   % gdad  % ������ ���� �����
KSI_x = del*atan2(T_x(1,2),T_x(2,2));               % gdad  % ������ ���� ��������
ELIR_angle_ = [TETA_x, FI_x, KSI_x];
end                                                 % End function ELIR_angle(q) %
% -----------------------------------------------------------------------
% ---------------- Begin ������� �������� ������������ � ���� ������ ----
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% ---------------- Begin ������� ������� �������� � �������������� ������������
% -----------------------------------------------------------------------
function T_qvant_ = T_qvant(q0, q1, q2, q3)         % Begin function T_qvant(q0, q1, q2, q3) %
% ������� ��������� ��������� ������� �������� � �������������� ������������
T_qvant_ = [(q0^2. + q1^2. - q2^2. - q3^2.)      (2.*(q1*q2 + q0*q3))          (2.*(q1*q3 - q0*q2));
                (2.*(q1*q2 - q0*q3))       (q0^2. - q1^2. + q2^2. - q3^2.)     (2.*(q2*q3 + q0*q1));
                (2.*(q1*q3 + q0*q2))           (2.*(q2*q3 - q0*q1))       (q0^2. - q1^2. - q2^2. + q3^2.);]; % �������
end                                                 % End function T_qvant(q0, q1, q2, q3) %
% -----------------------------------------------------------------------
% ---------------- End ������� ������� �������� � �������������� ������������
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% ---------------- Begin ������� ����������� �������������� �� ������� �
% -----------------------------------------------------------------------
function Lamda_ = Lamda_T(T_)                       % Begin function Lamda_T(T_) %
% ������� ��������� ����� ������������� �� �������
Lamda_0 = 0.5 * sqrt(1. + T_(1,1) + T_(2,2) + T_(3,3));
Lamda_1 = -(T_(3,2) - T_(2,3))/(4. * Lamda_0);
Lamda_2 = -(T_(1,3) - T_(3,1))/(4. * Lamda_0);
Lamda_3 = -(T_(2,1) - T_(1,2))/(4. * Lamda_0);
Lamda_ = [Lamda_0 Lamda_1 Lamda_2 Lamda_3];
end                                                 % End function Lamda_T(T_) %
% -----------------------------------------------------------------------
% ---------------- End ������� ����������� �������������� �� ������� �
% -----------------------------------------------------------------------

% 1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    8
% 2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  105
% 
% SGP4 TSINCE              X                Y                Z
% 
%       0.            2328.97048951   -5995.22076416    1719.97067261
%     360.00000000    2456.10705566   -6071.93853760    1222.89727783
%     720.00000000    2567.56195068   -6112.50384522     713.96397400
%    1080.00000000    2663.09078980   -6115.48229980     196.39640427
%    1440.00000000    2742.55133057   -6079.67144775    -326.38095856
% 
%                       XDOT             YDOT             ZDOT
% 
%                        2.91207230      -0.98341546      -7.09081703
%                        2.67938992      -0.44829041      -7.22879231
%                        2.44024599       0.09810869      -7.31995916
%                        2.19611958       0.65241995      -7.36282432
%                        1.94850229       1.21106251      -7.35619372
