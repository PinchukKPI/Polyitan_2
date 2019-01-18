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
    t.StartDelay = 3600; %60*60;
    t.Period = 3600;
    t.ExecutionMode = 'fixedDelay';
    t.TimerFcn = @(~,~)timer_done();
    start(t)
end

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
end




 