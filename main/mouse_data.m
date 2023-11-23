function mouse_data(post,sl,runstep1,runstep2)
%post must be either "_PC" of "_BF" to chose between
close all;

currdir=mfilename('fullpath');
coreDir=currdir(1:end-49);
dataDir = [coreDir,'\IBD_Mouse_MRIs\'];%
addpath(genpath(coreDir))


color=[0,0,1;1,0,0;0,1,0;0,0,0.172413793103448;1,0.103448275862069,0.724137931034483;1,0.827586206896552,0;0,0.344827586206897,0;0.517241379310345,0.517241379310345,1;0.620689655172414,0.310344827586207,0.275862068965517;0,1,0.758620689655172;0,0.517241379310345,0.586206896551724;0,0,0.482758620689655;0.586206896551724,0.827586206896552,0.310344827586207;0.965517241379310,0.620689655172414,0.862068965517241;0.827586206896552,0.0689655172413793,1;0.482758620689655,0.103448275862069,0.413793103448276;0.965517241379310,0.0689655172413793,0.379310344827586;1,0.758620689655172,0.517241379310345;0.137931034482759,0.137931034482759,0.0344827586206897;0.551724137931035,0.655172413793103,0.482758620689655;0.965517241379310,0.517241379310345,0.0344827586206897;0.517241379310345,0.448275862068966,0;0.448275862068966,0.965517241379310,1;0.620689655172414,0.758620689655172,1;0.448275862068966,0.379310344827586,0.482758620689655;0.620689655172414,0,0;0,0.310344827586207,1;0,0.275862068965517,0.586206896551724;0.827586206896552,1,0;0.724137931034483,0.310344827586207,0.827586206896552;0.241379310344828,0,0.103448275862069;0.931034482758621,1,0.689655172413793]
%SAMP
group = {'Stem','Cont','Dex','Norm'};
group_label={'Stem_S','Cont_S','Dex_S','Norm_S','Stem_T','Cont_T','Dex_T'}


%General name must pass in subtitled
Group_Stem= {'C822_rp','C822_np','C821_lp','C821_rp','1lp','1rp','np','2lp','2rp'};
Path_Stem = [14,10,16,14,8,7,5,8,9,nan,nan];
Group_Cont = {'C808_np','C822_llp','C823_lp','C823_rp','C823_np','5np','4lp','3lp','3rp','3np'};
Path_Cont = [10,16,12,18,18,10,7,6,3,5,nan];
Group_Dex = {'C821_lrp','C821_llp','C822_lrp','6lp','6rp','6np','7np'};
Path_Dex = [7,3,7,1,1,0,1,nan,nan,nan,nan];
Group_Norm = {'C141_lp','C141_np','C141_rp','C142_np'};
Path_Norm = [0,0,0,0,nan,nan,nan,nan,nan,nan,nan];

score=[Path_Stem;Path_Cont;Path_Dex;Path_Norm];
groupSamp=[1,1,1,1,0,0,0,0,0,0,0;1,1,1,1,1,0,0,0,0,0,0;1,1,1,0,0,0,0,0,0,0,0;1,1,1,1,0,0,0,0,0,0,0]';
tt(1,:)=[0,0,0,0,0,0];
groupSamp=groupSamp';
cd('D:Radiomics_Mouse_Analysis\IBD_Mouse_MRIs\Results');

featstack=[];
groupstack1=[];
groupstack2=[];
longestpatlist=11*3;%keep track of which sample has the most patients total
colorselect={};

P_Samp={};
feature_count=1;
which=[4,8,9,11]


if runstep1==1
for i=1:length(which)
    whichstat=which(i);
    for f=1:41
        mousecount=0;
        for g = 1:length(group);
            
            patients = eval(['Group_',group{g}]);
            count=1;
            feat=zeros(longestpatlist,1);
            groupSS=zeros(longestpatlist,1);
            for ii = 1:length(patients);
                mousecount=mousecount+1;

                load([dataDir,'\MAT',post,'\',patients{ii},post,'.mat']);
                
                temp1=(featstats(whichstat,f,:));
                temp1=reshape(temp1,[1,3]);
                feat(count+0)=temp1(1);
                feat(count+1)=temp1(2);
                feat(count+2)=temp1(3);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                groupSS(count+0)=groupSamp(g,ii);
                groupSS(count+1)=groupSamp(g,ii);
                groupSS(count+2)=groupSamp(g,ii);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                count=count+3;
                
                
                temp=eval(['Path_',group{g}]);
            end
            
            feat(feat==0) = NaN;
            feat=feat';
            
            groupSS=groupSS';
            
            
            featstack=[featstack;feat];
            groupstack1=[groupstack1;groupSS];
            clear feat groupSS groupTS
        end
        
        featCell{i,f}=featstack;
        group1Cell{i,f}=groupstack1;
 
        featstack=[];
        groupstack1=[];
    end
end
save(['MouseFeatStack',post,'.mat'],'featCell','group1Cell')
beep
else
    load(['MouseFeatStack',post,'.mat'])
    
end
if runstep2==1
for i=1:length(which)
        whichstat=which(i);
    for f=1:41
         mousecount=0;
        featstack=featCell{i,f};
        groupstack1=group1Cell{i,f};
        groupstack2=group2Cell{i,f};
        
        [P1,H1]=MouseTestingBank(featstack,groupstack1,1,.05,sl);
        [P2,H2]=MouseTestingBank(featstack,groupstack2,2,.05,sl);
        
        Samp_T=featstack(1,:).*groupstack1(1,:);
        Samp_T(Samp_T==0)=nan;
        Test_SvsT(:,1)=Samp_T;
        Test_SvsT(:,2)=TNF_T;
        PC=kruskalwallis(Test_SvsT(:,[1,2]),[],'off');
        PC=[whichstat,f,PC];
        P_Stem{i,f}=PC;
        
        Samp_T=featstack(2,:).*groupstack1(2,:);
        Samp_T(Samp_T==0)=nan;
        Test_SvsT(:,1)=Samp_T;
        PC=kruskalwallis(Test_SvsT(:,[1,2]),[],'off');
        PC=[whichstat,f,PC];
        P_Cont{i,f}=PC;
        
        Samp_T=featstack(3,:).*groupstack1(3,:);
        Samp_T(Samp_T==0)=nan;
        Test_SvsT(:,1)=Samp_T;
        PC=kruskalwallis(Test_SvsT(:),[],'off');
        PC=[whichstat,f,PC];
        P_Dex{i,f}=PC;
        
        P1=[whichstat,f,P1];
        P2=[whichstat,f,P2];
        H1=[whichstat,f,H1];
        H2=[whichstat,f,H2];
        
        P_Samp{i,f}=P1;
        H_Samp{i,f}=H1;
        
        feature_count=feature_count+1;
        
        featrr1=rescale_range(featstack,0,1)'; 
        groupstack1=groupstack1';
        
        featrr1(groupstack1==0)=nan;
        
        score1=score;
        score2=score;
        
        score1(groupSamp==0)=nan;
        score1=score1';
        
        if f==41
            a1=[num2str(f),'-',num2str(whichstat),'-int'];
            
        else
            a1=[num2str(f),'-','combined',num2str(4),'-',statnames{whichstat,f}];
            
        end
%         
% %         figure('Visible','Off','Position', [100, 100, 1024, 824])
% %         subplot(1,3,1);
% %         hold on
% %         boxplot(featrr1)
% %         title('Samp');
% %         ylim([0,1])
% %         set(gca,'XTickLabel',group_label(1:4));
% %         hold off
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         subplot(1,3,2);
% %         hold on
% %         boxplot(score1);
% %         title('Pathology score');
% %         ylim([0,20])
% %         set(gca,'XTickLabel',group_label(1:4));
% %         hold off
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %         subplot(1,3,3);
% %         hold on
% %         boxplot(score2);
% %         title('Pathology score');
% %         ylim([0,20])
% %         set(gca,'XTickLabel',group_label(5:7));
% %         hold off
% %         saveas(gcf,[a1,post,'_BP','.png'])
% %         close
%         
%         
%         xpos=ones(33,4);
%         xpos(:,2)=xpos(:,2)*2;
%         xpos(:,3)=xpos(:,3)*3;
%         xpos(:,4)=xpos(:,4)*4;
%         
%         xpos1=xpos;
%         xpos2=xpos;
%         
%         xpos1=reshape(xpos1,[1,33*4]);
%         xpos2=reshape(xpos2,[1,33*4]);
%         featrr1=reshape(featrr1,[1,33*4]);
%         featrr2=reshape(featrr2,[1,33*4]);
%         [featrr1,xpos1] =removeNaNs(featrr1,xpos1);
%         [featrr2,xpos2] =removeNaNs(featrr2,xpos2);
%         
%         
%         CreatePlot(color,a1,xpos1,featrr1,xpos2,featrr2,score1,score2,group_label,post)
%         clear feat featstack groupstack* groupTS groupSS

        
    end
end


P_Samp=P_Samp';
P_Stem=P_Stem';
P_Cont=P_Cont';
P_Dex=P_Dex';s
H_Samp=H_Samp';
save(['Feature_Performance',post,num2str(sl),'.mat'],'statnames','P_Stem','P_Samp','H_Samp','P_Cont','P_Dex','P_Samp')
end
