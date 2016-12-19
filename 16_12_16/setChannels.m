function [Chann] = setChannels
%% ask for the amount and channels names
% Chann=[channel names, Channel files]

aux=0;
while aux==0
    a1=input('number of channels? [insert a number] : ');
    
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
        %channel name (e.g.: GFP)
        Chann{i,1}=input(['name of the channel ',num2str(i),' : '],'s');
        
        %filename
        waitfor(msgbox(['assign a file for the channel ',num2str(i),',',char(10),'picking it in the next window'],'Warning'))
        Chann{i,2}=uigetfile({'*.xls;*xlsx','Excel Data(*.xls,*.xlsx)';'*.*',  'All Files (*.*)'},'Pick a file');

        %check
        if isequal(Chann{i,2},0)
            waitfor(msgbox(['User selected Cancel, please choose a file',char(10),'(re do all again)'], 'Error','error'))
        else
            disp(['User selected : ', Chann{i,2}])
            aux=1;
        end
    end
end
    
    
