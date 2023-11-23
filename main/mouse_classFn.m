clear
%code for training and testing 
addpath(genpath('D:\Radiomics_Mouse_Analysis\classification\main'))
addpath(genpath('D:\Radiomics_Mouse_Analysis\classifiers\IBD_Mouse_MRIs\Results'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('FeatSelCD.mat')
train_file=load('Sev2_Mouse_grp.mat')
test_file=load('Testing_Mouse_grp.mat')


train=train_file.MStack_grp;
test=test_file.MStack_grp;

train_class1=train{2};
train_class0=train{3};
label1=ones(size(train_class1,1),1);
label0=zeros(size(train_class0,1),1);

[train_class,mean_vec,mad_vec]=simplewhiten([train_class1;train_class0])

labels_full=[ones(1,24+20+16),zeros(1,14)]'

test_full=[test{1};test{2};test{3};test{4}];
test_full=simplewhiten(test_full,mean_vec,mad_vec);


sel=[spval_feat{1}(1,:);spval_feat{2}(1,:);spval_feat{3}(1,:);spval_feat{4}(1,:);spval_feat{5}(1,:);spval_feat{6}(1,:)]
sorted_sel=sortrows(sel,2)
sel=sorted_sel((1:3))

E3=TreeBagger(50,train_class(:,sel),[label1;label0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [label_full_rf3,score_full_rf3_12,cost_full_rf3]=predict(E3,test_full(:,sel));
 
set_score_rf3_12=floor(rescale_range(score_full_rf3_12(:,2),0,20))
score=score_full_rf3_12(:,2)
set=set_score_rf3_12;
save('RadScoresSAMP.mat','score','set','selfn','E3')
