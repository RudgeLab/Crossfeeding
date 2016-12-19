function [Labels] = setLabels
%% ask for the amount, channels names and data type
% Labels=[Label names, label files, data type]

aux=0;
while aux==0
    a1=input('number of labels? [insert a number] : ');
    
    %check if the input is valid
    if isnumeric(a1)==1
        if a1>=1
            aux=1;
        else waitfor(msgbox('invalid input, number have to be >=1','Warning','error'))
        end
    else waitfor(msgbox('invalid input, it have to be a number','Warning','error'))
    end
    %waitfor(msgbox()) == questdlg()
end
% the value is asigned only if the input is numeric and >=1.

for i=1:a1
    aux=0;
    while aux==0     %% only to check proper selection
        
        prompt = {'Label name','Data type (1=text, 2=num)'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'e.g: Strain / IPTG / etc','enter 1 or 2'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Labels{i,1}=answer{1};
        Labels{i,3}=str2num(answer{2});
        
        if Labels{i,3}==1
            disp(['label ',answer{1},' is text data'])
            aux=1;
        elseif Labels{i,3}==2
            disp(['label ',answer{1},' is numeric data'])
            aux=1;
        else
            waitfor(msgbox(['User put an invalid data type !',char(10),'(re do all again)'], 'Error','error'))
        end
    end
    aux2=0;
    while aux2==0
        %filename
        waitfor(msgbox(['assign a file for the label '' ',answer{1},' '' ,',char(10),'picking it in the next window'],'Warning'))
        Labels{i,2}=uigetfile({'*.xls;*xlsx','Excel Data(*.xls,*.xlsx)';'*.*',  'All Files (*.*)'},'Pick a file');
        
        %check
        if isequal(Labels{i,2},0)
            waitfor(msgbox(['User selected Cancel, please choose a file',char(10),'(select again)'], 'Error','error'))
        else
            disp(['User selected : ', Labels{i,2},char(10)])
            aux2=1;
        end
    end
end
end