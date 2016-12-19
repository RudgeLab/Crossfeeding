% This program read and analyse the data obtained from growing cells in
% 96 well plate reader.
% Specifiqly for ratiometric logic, it is, three channels are readed: CFP, YFP and OD

%Kevin things

function DataAnalysis(Ratio)%%[OD YFP CFP TM]=DataAnalysis
%Ratio=0 --> don´t perform ratiometric analysis
%Ratio = other value --> compute ratiometric

close all
Ratio=Ratio;
if isempty(Ratio)==1
    Ratio=input('do you want to compute ratiometric? (0 = No)');
end
%%%%%%%%%%%%%%%%%%
%Read and load Data

% Readed=dir('RawData.mat');
% if isempty(Readed)== 1   %check if the data is readed and saved already
%     ReadData            %read and save the data in RawData.mat file
% else
%     load RawData.mat  %load the saved data (with ReadData function)
% end
%%%%%
%read the files

Readed=dir('RawData.mat');

if isempty(Readed)== 0   %check if the data is readed and saved already
    
    load RawData.mat  %load the saved data (with ReadData function)
    
else
    CHs=setChannels();     %Call the function,  CHs=[Channel names, Channel files]
    CHnames =CHs(:,1);
    nfiles = size(CHs);
    
    %read the files
    for i=1:nfiles(1)
        [num(:,:,i),txt(:,:,i),raw(:,:,i)]=xlsread(CHs{i,2});
        % num = [time, well number, channels]
    end
    
    save('RawData.mat','num','txt','raw','CHnames')
    
end

%%%%

%% Define Labels
%%%%%

R_lab=dir('labels.mat');
if isempty(R_lab)== 0
    load labels.mat                 %If you wanna edit labels, open labels and edit on variable editor
else
    
    LabelsM=setLabels();               %call the function, LabelsM=[Label name, label filename, data type (1:text, 2: numeric)]
    nfileL = size(LabelsM);
    Lnames=LabelsM(:,1);
    Ltype=LabelsM(:,3);
    
    %read the files
    for i=1:nfileL(1)
        
        [Lnum(:,:,1),Ltxt(:,:,1),Lraw(:,:,1)]=xlsread(LabelsM{i,2});
        
        if LabelsM{i,3}==1              %1=text data
            C_labels{i}=strtrim(Ltxt(:,:,1));  %strtrim removes leading and trailing whitespace characters
        elseif LabelsM{i,3}==2          %2=numeric data
            C_labels{i}=Lnum(:,:,1);
        else
            waitfor(msgbox('The type of label data is not valid, check it please','Warning','error'))
            stop
        end
                
        %restart the variables (to avoid size matrix plomens)
        Lnum=[]; Ltxt=cell(0); Lraw=cell(0);
        
    end
    save('labels.mat','C_labels','Lnames','Ltype')
end
%%%%%%%%%%%%%%%%%%%%%%%
% Organize Data 
%%%%%

%get categories o replicates
[Rplic Groups LCat]=groupLab(C_labels); 
%LCat{2}
Rplic

nChan=size(num,3);      % number of channels
nLabels=length(Lnames); %number of labels

T=num(:,1,:)*24;           %Time vector for each channel (time data converted to hour units)
%%Thi is for splinetool
% nRep=8; %figure out how to know the number of replicates for each one
% for i=1:nRep
% TM(:,i)=T(:);  %Time Matrix, to match with each replicate (8) in splinetool
% end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%%% perform some preliminar plots

%% ask for the axis type
a1=input('use fixed axis? [Y/N] : ','s');
if strcmp(a1,'N')==1
    a1=0;
elseif strcmp(a1,'Y')==1
    a1=1;
else
    msgbox('invalid input, default value ''Y'' assigned','Warning')
    a1=1;
end
%%
%%%%


%% save max value of each dataset
for i=1:nChan
Ymax(i)=max(max(num(:,3:end,i)));  
end

if Ratio~=0
    %% cambiar este pedazo chanta por uno que pregunte al usuario por r1 y r2 (en base a la lista de los channels existentes)
    for i=1:size(CHnames,1)
        if strcmp(CHnames{i},'CFP')==1
            r2=i;
        end
        if strcmp(CHnames{i},'YFP')==1
            r1=i;
        end
    end
end
%%%%%
%% Background determination
%%%%%
% %plot OD background
% delay=11; % --> mediciones inciales en que condensación afecta lecturas
% startG=19;
% figure()
% plot(T(delay:startG,1),num(delay:startG,3:86,1),'o')
% hold on
% % ---> hay que tomar del 12 al 18 como blanco
OD_B=mean(mean(num(12:18,3:86,1)))
% plot([T(delay,1) T(startG,1)], [OD_B OD_B], 'r-')
delay=12;
% Blank substraction

num(:,3:end,1)=num(:,3:end,1)-OD_B;
%make > 0 if became negative
num(num<0)=0.0001;


%%
%perform some preliminar plots
for i=1:nChan
        figure(1)
        subplot(1,2,i)
        plot(T(delay:end,i),num(delay:end,3:end,i))
        title([CHnames{i}])
        xlabel('Time [hr]')
            
end

for i=1:nChan
    % la idea es que acá todos los graficos incluyen toda la data, pero en
    % cada label, cada catergoría tiene un color en particular.
        figure()
        for j=1:nLabels
        subplot(2,3,j)
        if j==1        
        end
        plot(T(delay:end,i),num(delay:end,3:end,i))
        title([Lnames{j}])
        end
        %legend([p1(1),p2(1),p3(1)],'SrpR','pAN1717','pAN1818','Location','northwest','FontSize',11)
        xlabel('Time [hr]')
end 



%plot threonine
figure()
hold on
cc=jet(7);
cch=hot(8);
ccw=winter(7);
T_IPTG=0;
for m=1:3
    %estos son los -IPTG
    if m==1 || m==2
        subplot(2,2,m)
        hold on
        for i=1:7
            Thre(i,:)=mn2m(i,7:2:12);
            p(i,:)=plot(T(delay:end,1),num(delay:end,Thre(i,:),1),'color',ccw(i,:));
        end
        Thre
        legend([p(1),p(2),p(3),p(4),p(5),p(6),p(7)],...
            [num2str(C_labels{1,4}(1,7)),' mM Thr'],[num2str(C_labels{1,4}(2,7)),' mM Thr'],...
            [num2str(C_labels{1,4}(3,7)),' mM Thr'],[num2str(C_labels{1,4}(4,7)),' mM Thr'],...
            [num2str(C_labels{1,4}(5,7)),' mM Thr'],[num2str(C_labels{1,4}(6,7)),' mM Thr'],...
            [num2str(C_labels{1,4}(7,7)),' mM Thr'],'Location','northwest','FontSize',11)
        if m==1
            title('-IPTG')
        end
    end
    %estos son los +IPTG
    if m==2 || m==3
        subplot(2,2,m)
        hold on
        for j=1:7
            Thre(j,:)=mn2m(j,8:2:12);
            p2(j,:)=plot(T(delay:end,1),num(delay:end,Thre(j,:),1),'color',cch(j,:));
        end
        Thre
        T_IPTG=Thre;
        legend([p2(1),p2(2),p2(3),p2(4),p2(5),p2(6),p2(7)],...
            [num2str(C_labels{1,4}(1,7)),' mM Thr'],[num2str(C_labels{1,4}(2,7)),' mM Thr'],...
            [num2str(C_labels{1,4}(3,7)),' mM Thr'],[num2str(C_labels{1,4}(4,7)),' mM Thr'],...
            [num2str(C_labels{1,4}(5,7)),' mM Thr'],[num2str(C_labels{1,4}(6,7)),' mM Thr'],...
            [num2str(C_labels{1,4}(7,7)),' mM Thr'],'Location','northwest','FontSize',11)
        if m==3
            title('+IPTG')
        else
            title('+/- IPTG')
        end
    end
end


figure()
hold on
T_IPTG(7,:) %solo toma la 7 para probar
%Td=T(delay:end,1);
plot(T(delay:end,1),num(delay:end,T_IPTG(7,:),1),'o')
%plot(T(delay:end,1),num(delay:end,T_IPTG(7,:),1),'color',cch(1:3,:))
Data=num(delay:end,T_IPTG(7,:),1);
Param=obtenerParam(Data,T(delay:end,1))

LP=size(Param,2);
for i=1:LP
y(:,i)=gompertz3(Param(:,i),T(delay:end,1));
ODt(:,i)=exp(y(:,i))*num(delay,T_IPTG(7,i),1);  %num(1,i) is the OD at first measure for every row
plot(T(delay:end,1),ODt(:,i))
end
legend('fit 1','fit 2','fit 3')

%%%%%%%%%%
%plot Isoleucina
figure()
hold on
I_IPTG=0;
for m=1:3
    %estos son los -IPTG
    if m==1 || m==2
        subplot(2,2,m)
        hold on
        for i=1:7
            ilv(i,:)=mn2m(i,1:2:6);
            p(i,:)=plot(T(delay:end,1),num(delay:end,ilv(i,:),1),'color',ccw(i,:));
        end
        ilv
        legend([p(1),p(2),p(3),p(4),p(5),p(6),p(7)],...
            [num2str(C_labels{1,3}(1,1)),' mM ilv'],[num2str(C_labels{1,3}(2,1)),' mM ilv'],...
            [num2str(C_labels{1,3}(3,1)),' mM ilv'],[num2str(C_labels{1,3}(4,1)),' mM ilv'],...
            [num2str(C_labels{1,3}(5,1)),' mM ilv'],[num2str(C_labels{1,3}(6,1)),' mM ilv'],...
            [num2str(C_labels{1,3}(7,1)),' mM ilv'],'Location','northwest','FontSize',11)
        if m==1
            title('-IPTG')
        end
    end
    %estos son los +IPTG
    if m==2 || m==3
        subplot(2,2,m)
        hold on
        for j=1:7
            ilv(j,:)=mn2m(j,2:2:6);
            p2(j,:)=plot(T(delay:end,1),num(delay:end,ilv(j,:),1),'color',cch(j,:));
        end
        ilv
        I_IPTG=ilv;
        legend([p2(1),p2(2),p2(3),p2(4),p2(5),p2(6),p2(7)],...
            [num2str(C_labels{1,3}(1,1)),' mM ilv'],[num2str(C_labels{1,3}(2,1)),' mM ilv'],...
            [num2str(C_labels{1,3}(3,1)),' mM ilv'],[num2str(C_labels{1,3}(4,1)),' mM ilv'],...
            [num2str(C_labels{1,3}(5,1)),' mM ilv'],[num2str(C_labels{1,3}(6,1)),' mM ilv'],...
            [num2str(C_labels{1,3}(7,1)),' mM ilv'],'Location','northwest','FontSize',11)
        if m==3
            title('+IPTG')
        else
            title('+/- IPTG')
        end
    end
end


% take maximum OD on +IPTG wells for ilv gradient
for i=1:7
%size(num(delay:end,ilv(i,:),1))
ilvMax(i,:)=max(num(delay:end,ilv(i,:),1))
end
%%%%%%%
maxIlvMean=mean(ilvMax(1:5,:),2);
maxIlvMean(6:7)=mean(ilvMax(6:7,2:3),2);
figure()
%plot(C_labels{1,3}(1:6,1),maxIlvMean(1:6),'r*',C_labels{1,3}(1:6,1),ilvMax(1:6,:),'bo')
plot(C_labels{1,3}(1:7,1),maxIlvMean,'r*',C_labels{1,3}(1:7,1),ilvMax,'bo')
legend('mean','well value')
xlabel('Isoleucine [mM]')
ylabel('OD_6_0_0')
title('maximum growth value per well (T7ilvCc)')
%%% Blanks

% Obtain TOP10 label value asociated
% for i=1:length(LCat{1})
%     if strcmp(LCat{1}(i),'TOP10')
%         nT10=i;
%     end
%     if strcmp(LCat{1}(i),'Blanco')
%         nB=i;
%     end
%     if strcmp(LCat{1}(i),'Blanco +Kan')
%         nBK=i;
%     end
% end
%%
% eventualmente hacer una funcion para esto de los blancos
% take fluorescente blank
c1=3;           %---> parte de 3 porque num parte de la 3era col
c2=1;
for i=1:size(Groups,1)
    for j=1:size(Groups,2)
        if (Groups(i,j,1)==nT10) == 1
            nNum(c2)=c1;
            c2=c2+1;
        end
        c1=c1+1;
    end
end

FBlank(:,:,:)=num(:,nNum,:); 
mFBlank=mean(FBlank,2);

%take OD blank
c1=3;
c2=1;
for i=1:size(Groups,1)
    for j=1:size(Groups,2)
        if (Groups(i,j,1)==nB) == 1 || (Groups(i,j,1)==nBK) == 1
            nNum2(c2)=c1;
            c2=c2+1;
        end
         c1=c1+1;
    end
end
ODBlank(:,:,:)=num(:,nNum2,:); 
mODBlank=mean(ODBlank,2);

%%
% Substract the blanks

%Do the matrix to coincide size in substraction
for i=1:size(num(:,3:end,:),2)
    MFBlank(:,i,:)=mFBlank;
    MODBlank(:,i,:)=mODBlank;
end

%OD
num(:,3:end,1)=num(:,3:end,1)-MODBlank(:,:,1);

%Fluor
num(:,3:end,2:3)=num(:,3:end,2:3)-MFBlank(:,:,2:3);

%correct negative values
num(num(:,:,:)<0)=0;

%%
% compute and plot the ratios
figure(2)
for j=3:size(num,2)
    R1=num(:,j,r1);
    R2=num(:,j,r2);     %R2
    plot(R2,R1)
    hold on
    
end
title([CHnames{r1},' / ', CHnames{r2}])
xlabel(CHnames{r2})
ylabel(CHnames{r1})

    for   j=3:size(num,2)
        Pdel=25;        % 25 ~ 4 hrs, 30 ~ 5hrs, no se toman los valores desde hora cero!
        R1=num(Pdel:end,j,r1);
        R2=num(Pdel:end,j,r2);
        
        P = polyfit(R2,R1,1);
        R(j-2)=P(1);              %slope save
        R0(j-2)=P(2);             %Coef. pos. save
        Ladj=polyval(P,R2);
        p3=plot(R2,Ladj,'k');
 %       legend([p1(1),p3(1)],'data','linear fit','Location','southeast')
    end
    
save('Ratio.mat','R','R0')

%% mean Ratios

for k=1:max(max(Rplic))
    c1=1;
    rep=1;
    for i=1:size(Rplic,1)
        for j=1:size(Rplic,2)
            if (Rplic(i,j)==k) == 1
                Rrep(k,rep)=R(c1);
                rep=rep+1;
            end
            c1=c1+1;
        end
    end
end

% Rplic
% for k=1:max(max(Rplic))
%     RG(:,k)=R(Rplic==k);
% end
% RG;
MRG=mean(Rrep,2);         %media de R entre replicas.
size(MRG);
%%
%plot the ratios

% eje X values
for i=1:length(LCat{2})
Xval(length(LCat{2})+1-i)=str2num(LCat{2}{i}); % lo inverti para que calzara
end
Xval(12)=0.1;       %asignado manual para grafico log

%     'LineWidth',2,...
%     'MarkerSize',10,

figure(3)
subplot(1,2,1)
% for i=1:12:(12*4)
% plot(Xval,MRG(i:(i+11)),'o')
% hold on
% end
plot(Xval,MRG(1:12),'b-',Xval,R(1:12),'bo',Xval,R(13:24),'bo','LineWidth',2,'MarkerSize',4)
hold on
plot(Xval,MRG(13:24),'r-',Xval,R(25:36),'ro',Xval,R(37:48),'ro','LineWidth',2,'MarkerSize',4)
plot(Xval,MRG(25:36),'g-',Xval,R(49:60),'go',Xval,R(61:72),'go','LineWidth',2,'MarkerSize',4)
%plot(Xval,R(73:84),'k-o')
title('Ratiometric','FontSize',15)    %TOp10
xlabel('IPTG [uM]','FontSize',12)
ylabel('Retiometric [no units]','FontSize',12)

%logplot
subplot(1,2,2)
p1=semilogx(Xval,MRG(1:12),'b-',Xval,R(1:12),'bo',Xval,R(13:24),'bo','LineWidth',2,'MarkerSize',4);
hold on
p2=semilogx(Xval,MRG(13:24),'r-',Xval,R(25:36),'ro',Xval,R(37:48),'ro','LineWidth',2,'MarkerSize',4);
p3=semilogx(Xval,MRG(25:36),'g-',Xval,R(49:60),'go',Xval,R(61:72),'go','LineWidth',2,'MarkerSize',4);
%p4=semilogx(Xval,R(73:84),'k-o')   %TOP10
title('Ratiometric','FontSize',15)
xlabel('IPTG [uM]','FontSize',12)
ylabel('Retiometric [no units]','FontSize',12)
%legend('SrpR','pAN1717','pAN1818','Top10')
legend([p1(1),p2(1),p3(1)],'SrpR','pAN1717','pAN1818','Location','northwest','FontSize',11)

% 
% for i=1:size(OD,3)
% figure(1)
% subplot(3,4,i)
% plot(T,OD(:,:,i))
% title(['Growth - ',C_labels{i}])
% xlabel('Time [hr]')
% ylabel('OD 600nm')
% if a1==1 ; ylim([0 Ymax(1)]); end
% 
% figure(2)
% subplot(3,4,i)
% plot(T,YFP(:,:,i))
% title(['YFP - ',C_labels{i}])
% xlabel('Time [hr]')
% ylabel('Intensity [AU]')
% if a1==1 ; ylim([0 Ymax(2)]); end
% 
% figure(3)
% subplot(3,4,i)
% plot(T,CFP(:,:,i))
% title(['CFP - ',C_labels{i}])
% xlabel('Time [hr]')
% ylabel('Intensity [AU]')
% if a1==1 ; ylim([0 Ymax(3)]); end
% 
% end

function [M]=mn2m(m,n)
%to convert m*n data into 96*1 data.
c=1;
for i=1:length(m)
    for j=1:length(n)
        M(c)=(m(i)-1)*12+n(j)+2;          % +2 because result data have 98 columns
        c=c+1;
    end
end

function Param = obtenerParam(M,T)
%Returns parameters fitted vector for every well
% M = Data Matrix, T= time vector
% Param = Parameters fitted with gompertz model

OD0=M(1,:);  %OD at first measurement point
ODt=M;  %matrix with OD at time t, for every well
L=size(OD0,2);

pIni=ones(3,L)*0.2;   % mu=0.2 it's take as start value for fiting, and it's idem to k and lambda

%%%%Computation of Log(OD(t)/OD0) from experimental data
for i=1:L
    
ODrt(:,i)=ODt(:,i)/OD0(i); %ODrt= OD(t)/OD0

end

LODrt=log(ODrt); %the model is fitted to this matrix


%%Fitting

%init variables
Poptim=zeros(3,L);
fval=zeros(1,L);
exitflag=zeros(1,L);  

for i=1:L
        %pIni(:,i)----> only to check!
[Poptim(:,i),fval(i),exitflag(i)] = fminsearch(@fo,pIni(:,i),[],T,LODrt(:,i));
    
end

Param=Poptim;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Minimal Square function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valor = fo(p,T,LODrt)

y=gompertz3(p,T);

valor =sum((y-LODrt).^2);
