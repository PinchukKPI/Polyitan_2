function star_time_ = sidereal_time(T, dUT);
% Результат в рвдианах 
%
% T [sec] - UTC time
% dUT [sec] - different between UTC and Sidereal Time
%
G = T(1);
M = T(2);
d = T(3);
hours = T(4);
min = T(5);
second = T(6);
%
JDres = JulianDay(T);
JD = JDres(1);
MJD = JDres(2);
UTCmjd = MJD + hours/24. + min/1440. + second/86400.;% - 0.125;
tc = UTCmjd + dUT/86400.;
Tu = (fix(tc) - 51544.5)/36525.;
Sm_o = 1.753368559233266 + (628.3319706888409 + (6.770714*10^(-6.) - 4.51*10^(-10.) * Tu)*Tu)*Tu;
r_s = 6.300388098984891 + (3.707456*10^(-10.) - 3.707*10^(-14.) * Tu)*Tu;
Sm = Sm_o + r_s*(tc - fix(tc));
%
% Дата в Юлианских столетиях Dj
% Юлианская дата начала эпохи (1900, январь 0, 12h)
JD0 = 2415020;
Dj = (JD - JD0 - 0.5)/36525;
%
% Звезное время на начало по Гринвичу текущих суток ts
ts_o = 1.7399368931356 + Dj*(628.331950990909 + (0.6755878064*10^(-5.))*Dj);
t = hours*60.*60. + min*60. + second + dUT;
ts = ts_o + 1.0027379093*t*2.*pi/86400;
%
% (ts - fix(ts/(2*pi))*2*pi)*180/pi
% (Sm - fix(Sm/(2*pi))*2*pi)
% (ts - fix(ts/(2*pi))*2*pi) - (Sm - fix(Sm/(2*pi))*2*pi)

star_time_ = ts;

% end
