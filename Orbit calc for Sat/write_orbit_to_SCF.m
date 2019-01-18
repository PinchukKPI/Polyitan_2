


PN_settings_T_orb_00_msec_UTC = uint8(0 );%0;
%     uint8_t PN_settings_T_orb_0_msec_UTC;//Время начальной точки орбиты в UTC (результат решения навигационной задачи)
PN_settings_T_orb_0_msec_UTC = uint8(0 );%0;

massive_to_write_int = [PN_settings_T_orb_00_msec_UTC PN_settings_T_orb_0_msec_UTC   ];
massive_to_write_int = typecast(swapbytes(massive_to_write_int) , 'uint16');
  
  %fprintf('%X \n',massive_to_write_int);
%     uint32_t PN_settings_T_orb_0_sec_UTC;//Время начальной точки орбиты в UTC (результат решения навигационной задачи)
PN_settings_T_orb_0_sec_UTC = uint32(1524124200);
%     float PN_settings_vx_0;              //X Составляющая Вектора скорости (результат решения навигационной задачи)
PN_settings_vx_0 = single(5096.97021484375);
%     float PN_settings_vy_0;              //Y Составляющая Вектора скорости (результат решения навигационной задачи)
PN_settings_vy_0 = single(3464.14404296875);
%     float PN_settings_vz_0;              //Z Составляющая Вектора скорости (результат решения навигационной задачи)
PN_settings_vz_0 = single(-4070.42700195313);
%     float PN_settings_x_0;               //X Координата (результат решения навигационной задачи)
PN_settings_x_0 = single(1095664.375);
%     float PN_settings_y_0;               //Y Координата (результат решения навигационной задачи)
PN_settings_y_0 = single(4320547.5);
%     float PN_settings_z_0;               //Z Координата (результат решения навигационной задачи)
PN_settings_z_0 = single(5054908.5);

PN_settings  = single([ typecast(PN_settings_T_orb_0_sec_UTC , 'single')  PN_settings_vx_0 PN_settings_vy_0 PN_settings_vz_0 PN_settings_x_0 PN_settings_y_0 PN_settings_z_0]);
massive_to_write_PN_settings = massive_to_write_single(PN_settings, 'uint8', 'uint16');
%fprintf('%X ',[massive_to_write_int massive_to_write_PN_settings]);
string_orbit_setting = num2str([massive_to_write_int massive_to_write_PN_settings], '%4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X %4.4X');
   

orbit_CSF = fopen('orbit_CSF.txt','w');
     
   
fprintf(orbit_CSF, '<?xml version="1.0" encoding="windows-1251"?>\r\n');
fprintf(orbit_CSF,'<Requests xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="Requests.xsd">\r\n');
fprintf(orbit_CSF,'<MB_req id="ID23" NameRec="1WorbitOpen" Description="открытие записи в настройки подсистемы навигации"><Mb_adr>1</Mb_adr><Mb_func>WRITE_16</Mb_func><Mb_pos>3787</Mb_pos><Mb_num>1</Mb_num><Mb_wr_data>A5A5</Mb_wr_data><Count>1</Count><Timeout_ms>100</Timeout_ms></MB_req><MB_req id="ID25" NameRec="1WorbitPN" Description="запись орбиты в настройки навигации"><Mb_adr>1</Mb_adr><Mb_func>WRITE_16</Mb_func><Mb_pos>3814</Mb_pos><Mb_num>15</Mb_num><Mb_wr_data>');
fprintf(orbit_CSF,string_orbit_setting);
fprintf(orbit_CSF,'</Mb_wr_data><Count>1</Count><Timeout_ms>100</Timeout_ms></MB_req><MB_req id="ID24" NameRec="1WorbitClose" Description="закрытие записи орбиты в настройки навигации"><Mb_adr>1</Mb_adr><Mb_func>WRITE_16</Mb_func><Mb_pos>3787</Mb_pos><Mb_num>1</Mb_num><Mb_wr_data>5A5A</Mb_wr_data><Count>1</Count><Timeout_ms>100</Timeout_ms></MB_req></Requests>\r\n');

fclose(orbit_CSF);