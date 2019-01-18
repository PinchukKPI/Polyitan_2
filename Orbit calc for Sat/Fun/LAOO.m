function LAOO_ = LAOO(Bc, Bo, Sc, So);
%
% ��� ������� ������� ������������ � ���� ��������.
% Bc - [��], ��� � ���
% Bo - [��], ��� � ���
% Sc - [-], ����������� �� ������ � ���
% So - [-], ����������� �� ������ � ���
%
Uc = Sc;

Vc_ = cross(cross(Sc,Bc),Sc);
Vc = Vc_/norm(Vc_);

Wc_ = cross(Sc,Bc);
Wc = Wc_/norm(Wc_);

Tc_ = [Uc,Vc,Wc];
%
Uo = So;

Vo_ = cross(cross(So,Bo),So);
Vo = Vo_/norm(Vo_);

Wo_ = cross(So,Bo);
Wo = Wo_/norm(Wo_);

To_ = [Uo,Vo,Wo];
%
Tco = Tc_*inv(To_);
%
LAOO_ = Tco;
%