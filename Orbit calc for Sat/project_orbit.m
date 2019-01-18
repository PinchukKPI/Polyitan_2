function varargout = project_orbit(varargin)
% PROJECT_ORBIT MATLAB code for project_orbit.fig
%      PROJECT_ORBIT, by itself, creates a new PROJECT_ORBIT or raises the existing
%      singleton*.
%
%      H = PROJECT_ORBIT returns the handle to a new PROJECT_ORBIT or the handle to
%      the existing singleton*.
%
%      PROJECT_ORBIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT_ORBIT.M with the given input arguments.
%
%      PROJECT_ORBIT('Property','Value',...) creates a new PROJECT_ORBIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project_orbit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project_orbit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project_orbit

% Last Modified by GUIDE v2.5 01-Jan-2015 00:58:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project_orbit_OpeningFcn, ...
                   'gui_OutputFcn',  @project_orbit_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before project_orbit is made visible.
function project_orbit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to project_orbit (see VARARGIN)

% Choose default command line output for project_orbit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

web_read_test()
% UIWAIT makes project_orbit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project_orbit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CALC_ORBIT_PUSHED.
function CALC_ORBIT_PUSHED_Callback(hObject, eventdata, handles)
% hObject    handle to CALC_ORBIT_PUSHED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = 'Зачекайте, йде розрахунок';
set(handles.text5,'String',str);
pause(0.2);
str1 = get(handles.first_string_NORAD, 'String');

str2 = get(handles.second_string_NORAD, 'String');


year = str2num(get(handles.year_text,'String'));

month = str2num(get(handles.month_text,'String'));

day = str2num(get(handles.day_text,'String'));

hours = str2num(get(handles.hours_text,'String'));

minutes = str2num(get(handles.minutes_text,'String'));

seconds = str2num(get(handles.sec_text,'String'));



str = orbit_two_lines_for_send_to_MS(str1, str2, year, month, day, hours, minutes, seconds);

set(handles.text5,'String',str);
%№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№№


function first_string_NORAD_Callback(hObject, eventdata, handles)
% hObject    handle to first_string_NORAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_string_NORAD as text
%        str2double(get(hObject,'String')) returns contents of first_string_NORAD as a double


% --- Executes during object creation, after setting all properties.
function first_string_NORAD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_string_NORAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function second_string_NORAD_Callback(hObject, eventdata, handles)
% hObject    handle to second_string_NORAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of second_string_NORAD as text
%        str2double(get(hObject,'String')) returns contents of second_string_NORAD as a double


% --- Executes during object creation, after setting all properties.
function second_string_NORAD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to second_string_NORAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function year_text_Callback(hObject, eventdata, handles)
% hObject    handle to year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of year_text as text
%        str2double(get(hObject,'String')) returns contents of year_text as a double


% --- Executes during object creation, after setting all properties.
function year_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function month_text_Callback(hObject, eventdata, handles)
% hObject    handle to month_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of month_text as text
%        str2double(get(hObject,'String')) returns contents of month_text as a double


% --- Executes during object creation, after setting all properties.
function month_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to month_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function day_text_Callback(hObject, eventdata, handles)
% hObject    handle to day_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of day_text as text
%        str2double(get(hObject,'String')) returns contents of day_text as a double


% --- Executes during object creation, after setting all properties.
function day_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to day_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hours_text_Callback(hObject, eventdata, handles)
% hObject    handle to hours_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hours_text as text
%        str2double(get(hObject,'String')) returns contents of hours_text as a double


% --- Executes during object creation, after setting all properties.
function hours_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hours_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minutes_text_Callback(hObject, eventdata, handles)
% hObject    handle to minutes_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minutes_text as text
%        str2double(get(hObject,'String')) returns contents of minutes_text as a double


% --- Executes during object creation, after setting all properties.
function minutes_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minutes_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sec_text_Callback(hObject, eventdata, handles)
% hObject    handle to sec_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sec_text as text
%        str2double(get(hObject,'String')) returns contents of sec_text as a double


% --- Executes during object creation, after setting all properties.
function sec_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sec_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on CALC_ORBIT_PUSHED and none of its controls.
function CALC_ORBIT_PUSHED_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to CALC_ORBIT_PUSHED (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
stop(t);
delete(t);
delete(hObject);

function web_read_test()
    global L1 L2 L1_arxiv L2_arxiv t
    % создаем файл архив на случай если его нет
    read_POLYITAN_arxiv = fopen('POLYITAN_arxiv.txt','a'); % Открытие файла *.txt
    fclose (read_POLYITAN_arxiv);
    % открываем архив на чтение
    read_POLYITAN_arxiv = fopen('POLYITAN_arxiv.txt','r'); % Открытие файла *.txt
    L1_arxiv = fgetl(read_POLYITAN_arxiv);
    L2_arxiv = fgetl(read_POLYITAN_arxiv);
    % вычитываем последние две строки
    while ~feof(read_POLYITAN_arxiv);
        L1_arxiv = fgetl(read_POLYITAN_arxiv);
        L2_arxiv = fgetl(read_POLYITAN_arxiv);
    end

    % смотрим на сайт
    read_url()

    % сравниваем вычитанную орбиту с архивом и сохраняем если нужно
    if ~strcmp(L1,L1_arxiv) || ~strcmp(L2,L2_arxiv)
        % создаем файл архив на случай если его нет
        read_POLYITAN_arxiv = fopen('POLYITAN_arxiv.txt','a'); % Открытие файла *.txt
        fprintf(read_POLYITAN_arxiv,'%s\r\n',L1);
        fprintf(read_POLYITAN_arxiv,'%s\r\n',L2);
        fclose(read_POLYITAN_arxiv);
    end

    L1_arxiv = L1;
    L2_arxiv = L2;

    % вычитываем каждый час сравниваем и если нужно дописываем в архив
    t = timer;
    t.StartDelay = 3600; %60*60; 3600
    t.Period = 3600; % 3600
    t.ExecutionMode = 'fixedDelay';
    t.TimerFcn = @(~,~)timer_done();
    start(t)

function timer_done()
    global L1 L2 L1_arxiv L2_arxiv
    disp('timer_done');
    read_url();
 
    % сравниваем вычитанную орбиту с архивом и сохраняем если нужно
    if ~strcmp(L1,L1_arxiv) || ~strcmp(L2,L2_arxiv)
        % создаем файл архив на случай если его нет
        read_POLYITAN_arxiv = fopen('POLYITAN_arxiv.txt','a'); % Открытие файла *.txt
        fprintf(read_POLYITAN_arxiv,'%s\r\n',L1);
        fprintf(read_POLYITAN_arxiv,'%s\r\n',L2);
        fclose(read_POLYITAN_arxiv);
        L1_arxiv = L1;
        L2_arxiv = L2;
    end
    
function read_url()
    global L1 L2
    % смотрим на сайт
    URL = 'http://www.celestrak.com/NORAD/elements/cubesat.txt';
    char_vector = urlread(URL,'Get',{'term','urlread'});
    %disp(char_vector);

    % создаем пустой файл
    cubesat = fopen('cubesat.txt','w');
    fprintf(cubesat,char_vector);
    fclose(cubesat);

    read_cubesat = fopen('cubesat.txt','r'); % Открытие файла *.txt

    string_read = fscanf(read_cubesat,'%s',1);
    % ищем наш спутник и две строки орбиты
    while strcmp(string_read,'') == 0
        if strcmp(string_read,'POLYITAN-2-SAU')

            thisline = fgetl(read_cubesat);
            L1 = fgetl(read_cubesat);
            L2 = fgetl(read_cubesat);
            %disp(L1);
            %disp(L2);
            break;
        end
        string_read = fscanf(read_cubesat,'%s',1);
    end
    fclose(cubesat);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global L1 L2
read_url();
set(handles.first_string_NORAD,'String',L1);
set(handles.second_string_NORAD,'String',L2);
 % handles.first_string_NORAD.String = L1;
 % handles.second_string_NORAD.String = L2;
 


