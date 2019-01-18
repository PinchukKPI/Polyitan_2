function GPStime_ = GPStime(UNIXtime);
%
% UNIXtime [-] - количество секунд с 1 января 1970
%
%
    GPS_UNIX_TIME_DIFF = 935280000;
    diff_GPS_UTC = 15;
    sumsec = (UNIXtime - GPS_UNIX_TIME_DIFF - diff_GPS_UTC);
    GPSweeks = fix(sumsec/604800);
    times = sumsec - GPSweeks * 604800;
    GPStime_ = [GPSweeks times];
end