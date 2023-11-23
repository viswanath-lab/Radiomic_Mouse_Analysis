clear%Split into training and testing data more clearly. 
post='_PC'
rrung1=0;
if rrung1==1

group = {'Stem','Cont','Dex','Norm'};
Group_Stem= {'822_1lp','822_rp','822_np','821_lp','821_rp','821_np'};
Group_Cont = {'808_np','822_llp','823_lp','823_lrp','823_np'};
Group_Dex = {'821_lrp','821_llp','822_lrp','828_np'};
Group_Norm = {'141_lp','141_np','141_rp','142_np'};

Path_Stem = [20,14,10,16,14,15];
Path_Cont = [10,16,12,18,18];
Path_Dex = [7,3,7,6];
Path_Norm = [0,0,0,0];

mices_nms={'822_1lp','822_rp','822_np','821_lp',...
    '821_rp','821_np','808_np','822_llp','823_lp',...
    '823_lrp','823_np','821_lrp','821_llp','822_lrp','828_np'...
    '141_lp','141_np','141_rp','142_np'}';

SelSlices_grp{1}=ones(6,4)
SelSlices_grp{2}=ones(5,4).*[NaN,NaN,NaN,NaN;1,1,1,1;NaN,NaN,NaN,NaN;1,1,1,1;1,1,1,1];
SelSlices_grp{3}=ones(4,4)
SelSlices_grp{4}=ones(4,4).*[1,1,1,1;NaN,1,1,NaN;1,1,1,NaN;1,1,1,1];%

which=[4,8,9,11];
trackmouse_grp={};
for g = 1:1:length(group);
        grp_stack=[];
    trackmouse={};
    count=1;
    patients = eval(['Group_',group{g}]);
    for ii = 1:length(patients);
        mouse_stack=[];
        load(['C:\Users\pvc5\Google Drive\Research\IBD_Mouse_MRIs\Mat_PC\',patients{ii},'-Axial',post,'.mat']);
        alpha=find(SelSlices_grp{g}(ii,:)==1);
        t_featstats=featstats(which,:,:);
        [a,b,c]=size(t_featstats);
        for ci=alpha
        trackmouse{end+1,1}=[patients{ii},num2str(ci)];
        end
        t_stack=[,t_featstats(1,:,alpha),t_featstats(2,:,alpha),t_featstats(3,:,alpha),t_featstats(4,:,alpha)];
        t_stack=reshape(t_stack,[b*4,length(alpha)])';
        mouse_stack=[mouse_stack,t_stack];
        grp_stack=[grp_stack;mouse_stack];
        trackmouse_grp{g}=trackmouse;
       
        
    end
    MStack_grp{g}=grp_stack;%This is to store the two groups
end
save('Sev2_Mouse_grp.mat','MStack_grp','trackmouse_grp')
load('Sev2_Mouse_grp.mat')
load('D:\Radiomics_Mouse_Analysis\classifiers\IBD_Mouse_MRIs\Mat_PC\822_1lp-Axial_PC.mat')
Max=MStack_grp{2};
Min=MStack_grp{3};
% stk=[Max;Min]
% [wnd_data,mean_vec,mad_vec] =simplewhiten(stk)
 [a,b]=size(Max);
% Max=wnd_data((1:a),:);
% Min=wnd_data((a+1):end,:);
pval=[];

for i=1:b
 p=ranksum(Max(:,i),Min(:,i))
 pval(i,1)=i;
 pval(i,2)=p;
end
sorted_p=sortrows(pval,2);
which=[4,8,9,11];
t_names=statnames(which,:);
names={}
for i=1:length(which)
    names=[names,t_names(i,:)]
end
sorted_names=names(sorted_p(:,1))'

for f=1:4
    figure
cell2boxplot(MStack_grp,sorted_p(f,1))
end

Greyfeat=[1:12,(1+102:12+102),(1+2*102:12+2*102),(1+3*102:12+3*102)]
Gradfeat=[13:25,(13+102:25+102),(13+2*102:25+2*102),(13+3*102:25+3*102)]
Harfeat=[26:64,(26+102:64+102),(26+2*102:64+2*102),(26+3*102:64+3*102)]
Colfeat=[89:101,(89+102:101+102),(89+2*102:101+2*102),(89+3*102:101+3*102)]
Gabfeat=[65:88,(65+102:88+102),(65+2*102:88+2*102),(65+3*102:88+3*102)]
Area=[102,204,306,408];

pval_grey=pval(Greyfeat,:)
pval_sbl=pval(Gradfeat,:)
pval_har=pval(Harfeat,:)
pval_col=pval(Colfeat,:)
pval_gab=pval(Gabfeat,:)
pval_area=pval(Area,:)

pval_feat{1}=pval_grey;
pval_feat{2}=pval_sbl;
pval_feat{3}=pval_har;
pval_feat{4}=pval_col;
pva_feat{5}=pval_gab;
pval_feat{6}=pval_area;


spval_grey=sortrows(pval(Greyfeat,:),2)
spval_sbl=sortrows(pval(Gradfeat,:),2)
spval_har=sortrows(pval(Harfeat,:),2)
spval_col=sortrows(pval(Colfeat,:),2)
spval_gab=sortrows(pval(Gabfeat,:),2)
spval_area=sortrows(pval(Area,:),2)

spval_feat{1}=spval_grey;
spval_feat{2}=spval_sbl;
spval_feat{3}=spval_har;
spval_feat{4}=spval_col;
spval_feat{5}=spval_gab;
spval_feat{6}=spval_area;

names_feat{1}=names(Greyfeat)
names_feat{2}=names(Gradfeat)
names_feat{3}=names(Harfeat)
names_feat{4}=names(Colfeat)
names_feat{5}=names(Gabfeat)
names_feat{6}=names(Area)


%%%%%%%%%%%%%%%%%%%%%%
%%%% full Comparison of all groups%%%%

c141_lp=[0,0,0,0];
c141_np=[1,0,0,0];
c141_rp=[0,0,0,0];
c142_np=[0,0,0,0];

c808_np=[2,2,1,1];
c821_llp=[1,1,0,0];
c821_lp=[2,2,1,2];
c821_lrp=[1,0,0,0];
c821_np=[2,1,2,0];
c821_rp=[1,2,1,0];
c822_1lp=[1,1,2,2];
c822_llp=[1,2,1,2];
c822_lrp=[1,1,0,0];
c822_np=[1,1,1,0];
c822_rp=[1,0,1,0];
c823_lp=[2,1,1,1];
c823_lrp=[1,1,2,0];
c823_np=[3,1,1,1];
c828_np=[1,0,1,1];

Samp_1=[0,0,0,0];
Samp_2=[0,0,0,0];
Samp_3=[2,2,2,2];
Samp_4=[1,1,2,1];
Samp_5=[3,2,2,1];
Samp_6=[1,1,2,2];
Samp_7=[0,0,1,1];
Samp_8=[0,1,0,0];
Samp_9=[0,0,1,0];
Samp_10=[0,1,0,0];
Samp_11=[2,1,0,1];

Rad_Stem= [c822_1lp';c822_rp';c822_np';c821_lp';c821_rp';c821_np'];
Rad_Cont = [c808_np';c822_llp';c823_lp';c823_lrp';c823_np'];
Rad_Dex = [c821_lrp';c821_llp';c822_lrp';c828_np'];
Rad_Norm=[c141_lp';c141_np';c141_rp';c142_np'].*[[1,1,1,1]';[NaN,1,1,1]';[NaN,1,1,1]';[1,1,1,1]']
Rad=[Rad_Stem;Rad_Cont;Rad_Dex;Rad_Norm]

save('FeatSelCD','pval','names','sorted_p','sorted_names','pval_feat','spval_feat','names_feat','Rad')
end
rrung4=0
if rrung4==1

%%%%%%%%%%%%%%%%%%%%%%%%
%                           SAMP3
Group = {'g3_SAMP1','g3_SAMP2','g3_SAMP3',...
    'g3_SAMP4','g3_SAMP5','g3_SAMP6',...
    'g3_SAMP7','g3_SAMP8','g3_SAMP9',...
    'g3_SAMP10','g3_SAMP11','g3_SAMP13',...
    'g3_SAMP14','g3_SAMP15','g3_SAMP16','g3_SAMP17'};

clear MStack_grp
which=[4,8,9,11];
trackmouse_grp={};
    grp_stack=[];
    trackmouse={};
    count=1;
    g=1
    patients = eval(['Group']);
       
    for ii = 1:length(patients);
        mouse_stack=[];
        load(['G:\My Drive\Research\IBD_Mouse_MRIs\Mat_PC\',patients{ii},'-Axial',post,'.mat']);
        
        t_featstats=featstats(which,:,:);
        [a,b,c]=size(t_featstats);
        for ci=1:c
        trackmouse{end+1,1}=[patients{ii},'-',num2str(ci)];
        trackmouse{end,2}=[featstats(4,102,ci)];
        end
         t_stack=[,t_featstats(1,:,:),t_featstats(2,:,:),t_featstats(3,:,:),t_featstats(4,:,:)];
        t_stack=reshape(t_stack,[b*length(which),c])';
        mouse_stack=[mouse_stack,t_stack];
        grp_stack=[grp_stack;mouse_stack];
        trackmouse_grp{g}=trackmouse;
        trackmouse_grp{1}=trackmouse;
        
        clear alpha
    end
    MStack_grp{1}=grp_stack;%This is to store the two groups
save('Testing_SAMP_Mouse_grp3.mat','MStack_grp','trackmouse_grp')

end
rrung5=1
if rrung5==1

%%%%%%%%%%%%%%%%%%%%%%%%
%                           SAMP4
Group = {'Samp-21-7','Samp-21-8','Samp-21-9','Samp-22-10','Samp-22-11','Samp-22-12','Samp-22-13','Samp-22-14','Samp-22-15','Samp-22-16','Samp-22-17','Samp-22-2','Samp-22-3','Samp-22-4','Samp-22-5','Samp-22-6','Samp-22-7','Samp-22-8','Samp-22-9'};

clear MStack_grp
which=[1,2,3,4];
trackmouse_grp={};
    grp_stack=[];
    trackmouse={};
    count=1;
    g=1
    patients = eval(['Group']);
       
    for ii = 1:length(patients);
        mouse_stack=[];
        load(['G:\My Drive\Research\IBD_Mouse_MRIs\Mat_PC\',patients{ii},'-Axial',post,'.mat']);
        
        t_featstats=featstats(which,:,:);
        [a,b,c]=size(t_featstats);
        for ci=1:c
        trackmouse{end+1,1}=[patients{ii},'-',num2str(ci)];
        trackmouse{end,2}=[featstats(4,102,ci)];
        end
         t_stack=[,t_featstats(1,:,:),t_featstats(2,:,:),t_featstats(3,:,:),t_featstats(4,:,:)];
        t_stack=reshape(t_stack,[b*length(which),c])';
        mouse_stack=[mouse_stack,t_stack];
        grp_stack=[grp_stack;mouse_stack];
        trackmouse_grp{g}=trackmouse;
        trackmouse_grp{1}=trackmouse;
        
        clear alpha
    end
    MStack_grp{1}=grp_stack;%This is to store the two groups
save('Testing_SAMP_Mouse_grp4.mat','MStack_grp','trackmouse_grp')

end
