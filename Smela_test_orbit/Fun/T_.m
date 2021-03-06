% -----------------------------------------------------------------------
% ---------------- Begin ������� ���������� ������� �������� ------------
% -----------------------------------------------------------------------
function T_x = T_(i,xx) % Begin function T_(i,xx) %
% ������� ��������� ��������� ������� �������� ������ �������� ��� ��
% �������� ����.
% �������� ������:
% i - ����� ���, ����� ��������� ��� ��������:
%       1 - ������� ������ ��� X,
%       2 - ������� ������ ��� Y,
%       3 - ������� ������ ��� Z;
% xx - �������� ���� �� ������� ���������� ���������, rad.
    if i == 1                                           % ������� ������ ��� X
        T_x = [ 1,        0,       0;
                0,  cos(xx), sin(xx);
                0, -sin(xx), cos(xx);];                 % ��������� �������
    end
    if i == 2                                           % ������� ������ ��� Y
        T_x = [cos(xx), 0, -sin(xx);
                     0, 1,        0;
               sin(xx), 0,  cos(xx);];                  % ��������� �������
    end
    if i == 3                                           % ������� ������ ��� Z
        T_x = [cos(xx), sin(xx), 0;
              -sin(xx), cos(xx), 0;
                     0,       0, 1;];                   % ��������� �������
    end

% end                     % End funktion T_(i,xx) %
% -----------------------------------------------------------------------
% ---------------- End ������� ���������� ������� �������� --------------
% -----------------------------------------------------------------------
