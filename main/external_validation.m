addpath(genpath('G:\My Drive\Research\IBD_Mouse_MRIs'))
load('RadScores_Samp_FNCD.mat')
load('FeatSelCDSAMP.mat')
sel=[spval_feat{1}(1,:);spval_feat{2}(1,:);spval_feat{3}(1,:);spval_feat{4}(1,:);spval_feat{5}(1,:);spval_feat{6}(1,:)]
sorted_sel=sortrows(sel,2)
sel=sorted_sel((1:3))


% %%%%%%%%%%%%%%%%%%%%%%%
% %group2
% test_file=load('Val_Mouse_grp.mat')
% test=test_file.MStack_grp;
% test_full=test;
% 
% 
% [label_full_rf3,score_full_rf3_grp2,cost_full_rf3]=predict(E3,test_full(:,sel));
% 
% score_grp2=score_full_rf3_grp2(:,2)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%
%group3
% test_file=load('Testing_SAMP_Mouse_grp3.mat')
% test=test_file.MStack_grp;
% test_full=test{1};
% 
% 
% [label_full_rf3,score_full_rf3_grp3,cost_full_rf3]=predict(E3,test_full(:,sel));
% 
% score_grp3=score_full_rf3_grp3(:,2)
% 
% save('SAMP_val.mat','score_grp2','score_grp3')
%%%%%%%%%%%%%%%%%%%%%%%
%group3
test_file=load('Testing_SAMP_Mouse_grp4.mat')
test=test_file.MStack_grp;
test_full=test{1};


[label_full_rf3,score_full_rf3_grp3,cost_full_rf3]=predict(E3,test_full(:,sel));

score_grp4=score_full_rf3_grp3(:,2)

save('SAMP_val.mat','score_grp4','-append')