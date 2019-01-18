function identification_angel_x = identification_angel(cos_x, sin_x)
    angel_asin = asin(sin_x);
    angel_acos = acos(cos_x);
    if cos_x > 0 & sin_x > 0
        % Первая четверть
        identification_angel_x = angel_asin;
    elseif cos_x < 0 & sin_x > 0
        % вторая четверть
          identification_angel_x = angel_acos;
    elseif cos_x < 0 & sin_x < 0
        % третья четверть
        identification_angel_x = 2.*pi - angel_acos;
    elseif cos_x > 0 & sin_x < 0
        % четвертая четверть
        identification_angel_x = 2.*pi - angel_acos; %SO 2.02.2013: fix the bug %angel_asin;
    elseif angel_asin == angel_acos
        identification_angel_x = angel_acos;
    end 
% end

