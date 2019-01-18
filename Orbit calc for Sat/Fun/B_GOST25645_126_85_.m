% -----------------------------------------------------------------------
% ---------------- Begin ������� ������� ������� �������� ��� -----------
% -----------------------------------------------------------------------
function B_GOST25645_126_85 = B_GOST25645_126_85_(F1, L, H1, DateTime, YEAR_GOST)           % Begin function B_GOST25645_126_85(x, y, z) %
% ������� ��������� ��������� ������ �������� ��� �� ��������� ������
% �������� ������:
% F1 - ������������� (�������������) ������ ����� � ������������ [���];
% L - ������������� (�������������) ������� ����� � ������������ [���];
% H1 - ������ ����� ��� ������� ���� [��];
% DateTime = [Year Month Day Hour Minute Second] - ������� ���� � ����� (������������� ������: Year Month Day);
% YEAR_GOST - ���� ����������� ������������� �������������;
%
% �������� ��������:
% X - �������� ������� �������� �� ��� x, ������������ �� ��������������� ��������� (�� �����) [���];
% Y - �������� ������� �������� �� ��� y, ������������ �� ��������� (�� ������) [���];
% Z - �������� ������� �������� �� ��� z, ������������ ����������� ���� [���];

% clc
% format long

% % �������1 1989
% H1 = 100
% L = 58*pi/180
% F1 = 80.6*pi/180
% % �����: 
% X_GOST = 4507.0;
% Y_GOST = 2446.9;
% Z_GOST = 53865.8;
% 
% % �������2 1989
% H1 = 3000
% L = 58*pi/180
% F1 = 80.6*pi/180
% % �����: 
% X_GOST = 2109.1;
% Y_GOST = -127.9;
% Z_GOST = 18644.0;
% 
% �������3 1989
% H1 = 6371
% L = 58*pi/180
% F1 = 80.6*pi/180
% % �����: 
% X_GOST = 940.1;
% Y_GOST = -191.9;
% Z_GOST = 7462.1;
% 
% �������4 1989
% H1 = 6385
% L = 58*pi/180
% F1 = 80.6*pi/180
% % �����: 
% X_GOST = 937.3;
% Y_GOST = -191.6;
% Z_GOST = 7437.9;
% 
% �������5 1989
% H1 = 12742.4
% L = 58*pi/180
% F1 = 80.6*pi/180
% % �����: 
% X_GOST = 298.4;
% Y_GOST = -90.8;
% Z_GOST = 2199.9;
% 
% �������6 1989
% H1 = 40000
% L = 58*pi/180
% F1 = 80.6*pi/180
% % �����: 
% X_GOST = 21.8;
% Y_GOST = -9.2;
% Z_GOST = 151.7;
% 
%     KOR_TO_GEOD = KOR_TO_GEOD_(x, y, z);
%     F1 = KOR_TO_GEOD(1); % B
%     L = KOR_TO_GEOD(2); % L
%     H1 = KOR_TO_GEOD(3); % H

    % ��������� �� ���� 25645.126-85
%     YEAR = 1989; % �������� ����
    if YEAR_GOST == 1985
        NH = 6; % ����� ��������
        K = (NH*NH + 3*NH)/2; % ���������� �������������
        G = [-29877, -1903, -2073, 3045, 1691, 1300, -2208, 1244, 835, 937, 780, 363, -426, 169, -215, 356, 253, -94, -161, -48, 52, 65, 50, -186, 4, 17, -102, 75, -61, 2, 24, -6, 4, 9, 0, 21, 6, 0, -11, -9, 2, 4, 4, -6, 5, 10, 1, -12, 9, -3, -1, 7, 2, -5, -4, -4, 2, -5, -2, 5, 3, 1, 2, 3, 0];
        G1 = [0, 5497, 0, -2191, -309, 0, -312, 284, -296, 0, 233, -250, 68, -298, 0, 47, 148, -155, -75, 95, 0, -16, 90, 69, -50, -4, 20, 0, -82, -26, -1, 23, 17, -21, -6, 0, 7, -21, 5, -25, 11, 12, -16, -10, 0, -21, 16, 9, -5, -6, 9, 10, -6, 2, 0, 1, 0, 3, 6, -4, 0, -1, 4, 0, -6];
        DG = [19.7, 11.5, -12.6, 1.8, 1.4, 4.3, -6.1, -0.7, -3.8, -0.4, 0.2, -7.4, -0.4, -5.7, 1.2, -0.1, -1.2, -2.4, -0.3, 0.5, 1.4, -0.4, 1.6, 0.9, -0.1, 0.7, 1.0, 0.4, -0.6, -0.1, 0.2, 0.9, 0.9, 0.3, 1.0, 0.5, -0.3, -0.1, 0.6, -0.7, 0.1, 0.2, -0.9, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        DG1 = [0, -20, 0, -16.4, -15.9, 0, 4.6, 2.8, -9.8, 0, 3.5, 2, 3.7, -0.3, 0, 0, 0.7, 0.1, 1.1, -0.1, 0, -0.7, -1.2, -0.3, -1.3, 0.4, 1.1, 0, 1, 0.3, 0.8, 0.7, -0.2, 0.3, 0, 0, 0.6, -0.3, 0.4, 0, 0.6, -1.2, -0.1, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    elseif YEAR_GOST == 2010
        NH = 6; % ����� ��������
        K = (NH*NH + 3*NH)/2; % ���������� �������������
        G = [-29496.5, -1585.9, -2396.6, 3026.0, 1668.6, 1339.7, -2326.3, 1231.7, 634.2, 912.6, 809.0, 166.6, -357.1, 89.7, -231.1, 357.2, 200.3, -141.2, -163.1, -7.7, 72.8, 68.6, 76.0, -141.4, -22.9, 13.1, -77.9, 80.4, -75.0, -4.7, 45.3, 14.0, 10.4, 1.6, 4.9, 24.3, 8.2, -14.5, -5.7, -19.3, 11.6, 10.9, -14.1, -3.7, 5.4, 9.4, 3.4, -5.3, 3.1, -12.4, -0.8, 8.4, -8.4, -10.1, -2.0, -6.3, 0.9, -1.1, -0.2, 2.5, -0.3, 2.2, 3.1, -1.0, -2.8, 3.0, -1.5, -2.1, 1.6, -0.5, 0.5, -0.8, 0.4, 1.8, 0.2, 0.8, 3.8, -2.1, -0.2, 0.3, 1.0, -0.7, 0.9, -0.1, 0.5, -0.4, -0.4, 0.2, -0.8, 0.0, -0.2, -0.9, 0.3, 0.4, -0.4, 1.1, -0.3, 0.8, -0.2, 0.4, 0.0, 0.4, -0.3, -0.3];
        G1 = [0, 4945.1, 0, -2707.7, -575.4, 0, -160.5, 251.7, -536.8, 0, 286.4, -211.2, 164.4, -309.2, 0, 44.7, 188.9, -118.1, 0.1, 100.9, 0, -20.8, 44.2, 61.5, -66.3, 3.1, 54.9, 0, -57.8, -21.2, 6.6, 24.9, 7.0, -27.7, -3.4, 0, 10.9, -20.0, 11.9, -17.4, 16.7, 7.1, -10.8, 1.7, 0, -20.5, 11.6, 12.8, -7.2, -7.4, 8.0, 2.2, -6.1, 7.0, 0, 2.8, -0.1, 4.7, 4.4, -7.2, -1.0, -4.0, -2.0, -2.0, -8.3, 0, 0.1, 1.7, -0.6, -1.8, 0.9, -0.4, -2.5, -1.3, -2.1, -1.9, -1.8, 0, -0.8, 0.3, 2.2, -2.5, 0.5, 0.6, 0.0, 0.1, 0.3, -0.9, -0.2, 0.8, 0, -0.8, 0.3, 1.7, -0.6, -1.2, -0.1, 0.5, 0.1, 0.5, 0.4, -0.2, -0.5, -0.8];
        DG = [11.4, 16.7, -11.3, -3.9, 2.7, 1.3, -3.9, -2.9, -8.1, -1.4, 2.0, -8.9, 4.4, -2.3, -0.5, 0.5, -1.5, -0.7, 1.3, 1.4, -0.3, -0.3, -0.3, 1.9, -1.6, -0.2, 1.8, 0.2, -0.1, -0.6, 1.4, 0.3, 0.1, -0.8, 0.4, -0.1, 0.1, -0.5, 0.3, -0.3, 0.3, 0.2, -0.5, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        DG1 = [0, -28.8, 0, -23.0, -12.9, 0, 8.6, -2.9, -2.1, 0, 0.4, 3.2, 3.6, -0.8, 0, 0.5, 1.5, 0.9, 3.7, -0.6, 0, -0.1, -2.1, -0.4, -0.5, 0.8, 0.5, 0, 0.6, 0.3, -0.2, -0.1, -0.8, -0.3, 0.2, 0, 0.0, 0.2, 0.5, 0.4, 0.1, -0.1, 0.4, 0.4, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    end
    YEAR1 = (DateTime(3)/30 +DateTime(2)-1)/12 + DateTime(1) - YEAR_GOST;
    for I = 1:K
        G(I) = G(I) + DG(I)*YEAR1;
        G1(I) = G1(I) + DG1(I)*YEAR1;
    end
    
    RS = 6371.200; % ������� ������ �����, ��
    A = 6378.200; % ������� ������� ������� ���������� ��������, ��
    A2 = A*A;
    A4 = A2*A2;
    B = 6356.800; % ����� ������� ������� ���������� ��������, ��
    B2 = B*B;
    B4 = B2*B2;
    A3 = 0.0000001;
    NH = NH + 1;
    
    
    cos2FI = (cos(F1))^2;
    sin2FI = (sin(F1))^2;
    sum1 = A2*cos2FI + B2*sin2FI;
    sum2 = A4*cos2FI + B4*sin2FI;
    R1 = (H1*H1 + 2*H1*(sum1)^0.5 + sum2/sum1)^0.5;
    F3 = atan(tan(F1)*(B2 + H1*sum1^0.5)/(A2 + H1*sum1^0.5));
    
    F = 0.5*pi - F3;
    S1=sin(F1-F3);
    S2=cos(F1-F3);
    
    C1=sin(F);
    C2=cos(F);
    S(1,1) = 1;
    
    for N = 2:NH
        S(1,N) = S(1,N - 1)*(2*N - 3)/(N - 1);
        S(2,N) = S(1,N)*sqrt((N - 1)*2/N);
        if N >= 3
            for M = 3:N
                S(M,N)= S(M-1,N)*sqrt((N-M+1)/(N+M-2));
            end
        end
    end
    
    P(1,1) = 1; P(1,2) = C2; P(2,2) = C1;
    R(1,1) = 0; R(1,2) = -C1; R(2,2) = C2;
    for N = 3:NH
        for M = 1:N
            if (M-N) < 0
                NR = ((N - 2)^2 - (M - 1)^2)/((2*N - 3)*(2*N - 5));
                P(M,N) = C2*P(M,N - 1) - NR*P(M,N - 2);
                R(M,N) = C2*R(M,N - 1) - C1*P(M,N - 1) - NR*R(M,N - 2);
            elseif (M-N) == 0
                P(M,N) = C1*P(M - 1,N - 1);
                R(M,N) = C1*R(M - 1,N - 1) + C2*P(M - 1,N - 1);
            else
                P(M,N) = 0;
                R(M,N) = 0;
            end
        end
    end
    for N = 1:NH
        for M = 1:N
            P(M,N) = P(M,N)*S(M,N);
            R(M,N) = R(M,N)*S(M,N);
        end
    end
    
    
    for M = 1:NH
        U1(M) = sin((M-1)*L);
        U2(M) = cos((M-1)*L);
    end      
    L1 = RS/R1;     
    A1 = abs(sin(F));
    if A1 < A3
        A1 = A3;
    else
        A1 = sin(F);
    end
    X = 0;
    Y = 0;
    Z = 0;
    J = 0;
    for N = 2:NH
        for M = 1:N
            A2 = (M - 1)/A1;
            J = J + 1;
            X = X +        (G(J)*U2(M) + G1(J)*U1(M))*L1^(N + 1)*R(M,N);
            Y = Y +        (G(J)*U1(M) - G1(J)*U2(M))*L1^(N + 1)*P(M,N)*A2;
            Z = Z + (-1)*N*(G(J)*U2(M) + G1(J)*U1(M))*L1^(N + 1)*P(M,N);
        end
    end
    %fprintf('G[6]   %12.8f   %12.8f   %12.8f  %12.8f   %12.8f   %12.8f   \n',G(1), G(2), G(3), G(4), G(5), G(6) );
    %fprintf('G1[6]  %12.8f   %12.8f   %12.8f  %12.8f   %12.8f   %12.8f   \n',G1(1), G1(2), G1(3), G1(4), G1(5), G1(6) );
    %X
    %Y
    %Z
    
    X = X*S2+Z*S1;
    Z = Z*S2-X*S1;
    
    B_GOST25645_126_85 = [X Y Z];
%     
%     [X Y Z]
%     [100*(X-X_GOST)/X_GOST 100*(Y-Y_GOST)/Y_GOST 100*(Z-Z_GOST)/Z_GOST]
%
end                                         % End function B_GOST25645_126_85(x, y, z) %
% -----------------------------------------------------------------------
% ---------------- End ������� ������� ������� �������� ��� -------------
% -----------------------------------------------------------------------
 