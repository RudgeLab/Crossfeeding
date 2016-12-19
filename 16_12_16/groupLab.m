function [GG,Group,L]=groupLab(Labels)
%% 
%GG = all the replicates groups (i,e: same number implies replicates)
%GG = [(A,B..H),(1...12)]
%Group = same number implies same group inside the label, 
%Group =[(A,B..H),(1...12),label]
%L=list of categories inside each label , L=[categories, Label] ->cellArray

S=length(Labels);
S1=size(Labels{1},1);
S2=size(Labels{1},2);

for i=1:S
    %found the categories
    CatL{i}=categorical(Labels{i});
    Laux=categories(CatL{i});
    L{i}=Laux;          %save it on L{i}
    %% this section is to know the group of each one element in label
    if isnumeric(Labels{i})
        for h=1:length(Laux)
            Laux{h}=str2num(Laux{h}); 
        end
    end
    % the above if, is to correct the data format (if is numeric, then, put
    % they as numbers)
    
    for k=1:length(Laux)
        for m=1:S1
            for n=1:S2
                if isnumeric(Labels{i})
                    if (Laux{k}==Labels{i}(m,n))==1
                        Group(m,n,i)=k;
                    end
                else
                    if strcmp(Laux{k},Labels{i}{m,n})==1
                        Group(m,n,i)=k;
                    end
                end
            end
        end
    end
    %%
end
%Group(:,:,:)

c=1;
%GZ=zeros(S1,S2);    %auxiliar to compare
GG=zeros(S1,S2);    % is the final grouped group
for i=1:S1
    for j=1:S2
        
        if i==1 && j==1
            G(c,:)=Group(i,j,:);
            %G2(
            n=1; m=1;
            for m=1:S1
                for n=1:S2
                    Group_aux(1,:)=Group(m,n,:);
                    if prod(G(c,:)==Group_aux)==1
                        GG(m,n)=c;
                    end
                end
            end
        else
%                         i
%                         j
            Group_aux(1,:)=Group(i,j,:);
            %             if i==1 && j==3
            %                 G(c,:)==Group_aux(1,:)
            %             end
            existG=0;
            comAll=[];
            for h=1:size(G,1)
                comAll(:,h)=(G(h,:)==Group_aux(1,:));
            end
            pg=prod(comAll,1);
            for h=1:size(pg,2)
                if pg(h) == 1
                    existG=1;
                end
            end
            if existG == 0      %only if is not a group recorded before
                %            if prod(G(c,:)~=Group_aux(1,:)) == 0        
                c=c+1;
                G(c,:)=Group(i,j,:);
                n=1; m=1;
                for m=1:S1
                    for n=1:S2
                        Group_aux(1,:)=Group(m,n,:);
                        if prod(G(c,:)==Group_aux)==1
                            GG(m,n)=c;

                        end
                    end
                end
            end
        end
        
        
    end

end
%size(GG)
%GG

