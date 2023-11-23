function Cn=CreatePlot(color,a1,xpos1,featrr1,xpos2,featrr2,score1,score2,group_label,post)
%load('Feature_Performance','Group_Color')
figure('Visible','Off','Position', [100, 100, 1024, 824])
subplot(2,2,1);
hold on
color1=color((1:16),:);
count=0;
GroupLeg={};
for i=1:length(featrr1)
    if mod(i,3)==1
        count=count+1;
    end
    if(isnan(featrr1))
    else
        scatter(xpos1(i),featrr1(i),74,color1(count,:),'filled')
    end
end
title('Samp');
xlim([0,5])
ylim([0,1])
%legend(GroupLeg)
set(gca,'XTick',(1:1:4),'XTickLabel',group_label(1:4));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2);
hold on
boxplot(score1);
title('Pathology score');
ylim([0,20])
set(gca,'XTickLabel',group_label(1:4));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
hold on
color2=color((16:end),:);
count=0;
for i=1:length(featrr2)
    if mod(i,3)==1
        count=count+1;
    end
    
    if(isnan(featrr1))
    else
        scatter(xpos2(i),featrr2(i),74,color2(count,:),'filled')
    end
end
title('TNF');
ylim([0,1])
xlim([0,4])
set(gca,'XTick',(1:1:3),'XTickLabel',group_label(5:7));
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4);
hold on
boxplot(score2);
title('Pathology score');
ylim([0,20])
set(gca,'XTickLabel',group_label(5:7));
hold off
saveas(gcf,[a1,post,'.png'])
close

end