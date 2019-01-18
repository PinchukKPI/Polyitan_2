function JulianDay_ = JulianDay(T);

%Input data
G = T(1);
M = T(2);
d = T(3);
hours = T(4);
min = T(5);
second = T(6);
%
% Значение Юлианской даты JD
JD = d + fix(1461.*(G + 4800. + fix((M - 14)/12.))/4. + fix(367.*(M - 2. - fix((M - 14.)/12.)*12.)/12.) - fix(3.*fix((G + 4900. + fix((M - 14.)/12.))/100.)/4.) - 32075);
%
year=G-1900;
month=M-3;
if month<0
    month=month+12.;
    year=year-1;
end
JDmjd = 15078.+365.*year+fix(year*0.25)+fix(0.5+30.6*month)+d;
%
MJD = JDmjd;% - 2400000.5;
%
JulianDay_ = [JD MJD];
%
