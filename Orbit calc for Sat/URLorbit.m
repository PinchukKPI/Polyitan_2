function varargout = URLorbit(varargin)
% URLORBIT MATLAB code for URLorbit.fig
%      URLORBIT, by itself, creates a new URLORBIT or raises the existing
%      singleton*.
%
%      H = URLORBIT returns the handle to a new URLORBIT or the handle to
%      the existing singleton*.
%
%      URLORBIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in URLORBIT.M with the given input arguments.
%
%      URLORBIT('Property','Value',...) creates a new URLORBIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before URLorbit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to URLorbit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help URLorbit

% Last Modified by GUIDE v2.5 01-Jan-2015 02:22:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @URLorbit_OpeningFcn, ...
                   'gui_OutputFcn',  @URLorbit_OutputFcn, ...
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


% --- Executes just before URLorbit is made visible.
function URLorbit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to URLorbit (see VARARGIN)

% Choose default command line output for URLorbit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes URLorbit wait for user response (see UIRESUME)
% uiwait(handles.figure1);

web_read_test()

% --- Outputs from this function are returned to the command line.
function varargout = URLorbit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global t
stop(t)
delete(t);


%str = 'Closing';
%set(handles.text1,'String',str);
%pause(2);

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

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





 
