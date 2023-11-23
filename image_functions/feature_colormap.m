function [] = feature_colormap(img1_path,img2_path,savedir,fselect,sufix)%mouse
cd(savedir)
a='low'
b='high'
%fselect=[80,30]
load(img1_path);
feat_Vol1=features;
Vol1=imgvol;
Mask1=maskvol;
Mask1(Mask1==1)=0;
Mask1(Mask1==2)=1;
Mask1(Mask1==4)=0;
S1=slices;
clear features imgvol maskvol slices

load(img2_path);
feat_Vol2=features;
Vol2=imgvol;
Mask2=maskvol;

Mask2(Mask2==1)=0;
Mask2(Mask2==2)=1;
Mask2(Mask2==4)=0;%check mask to make sure only a single mask
S2=slices;


setslices=3%vary this to slice you want to view
for j=1:length(fselect)
    if setslices>0
        for i=1:setslices

            img1=[Vol1(:,:,S1(i))];
            img2=[Vol2(:,:,S2(i))];
            % img3=[Vol3(:,:,S3(i))];

            img1rr=rescale_range(img1,0,1);
            img2rr=rescale_range(img2,0,1);
            %img3rr=rescale_range(img3,0,1);

            gt1=graythresh(img1rr);
            gt2=graythresh(img2rr);
            %gt3=graythresh(img3rr);

            tot_mask1=imbinarize(img1rr,gt1);
            tot_mask2=imbinarize(img2rr,gt2);
            % tot_mask3=imbinarize(img3rr,gt3);

            % tot_mask1=imfill(tot_mask1,'holes');
            %tot_mask2=imfill(tot_mask2,'holes');
            % tot_mask3=imfill(tot_mask3,'holes');

            m1=[Mask1(:,:,S1(i))];
            m2=[Mask2(:,:,S2(i))];
            % m3=[Mask3(:,:,S3(i))];

            m1(m1==1)=0;
            m2(m2==1)=0;
            % m3(m3==1)=0;
            m1(m1==3)=0;
            m2(m2==3)=0;
            %  m3(m3==3)=0;
            m1(m1==4)=0;
            m2(m2==4)=0;
            %  m3(m3==4)=0;
            m1(m1==2)=1;
            m2(m2==2)=1;
            % m3(m3==2)=1;

            cmap=colormap('jet');
            coloraxis=[0,1]

            img1(tot_mask1==0)=nan;
            img2(tot_mask2==0)=nan;
            %img3(tot_mask3==0)=nan;


            stack_img(:,:,i)=img1;
            stack_img(:,:,i+4)=img2;
            %        stack_img(:,:,i+8)=img3;

            f1=[feat_Vol1(:,:,i,fselect(j))];
            f2=[feat_Vol2(:,:,i,fselect(j))];
            %f3=[feat_Vol3(:,:,i,fselect(j))];


            f1(tot_mask1==0)=nan;
            f2(tot_mask2==0)=nan;
            %   f3(tot_mask3==0)=nan;

            stack_feat(:,:,i)=f1;
            stack_feat(:,:,i+4)=f2;
            %  stack_feat(:,:,i+8)=f3;

            stack_mask(:,:,i)=m1;
            stack_mask(:,:,i+4)=m2;
            %   stack_mask(:,:,i+8)=m3;

        end
    else
        stack_img=cat(3,Vol1(:,:,S1),Vol2(:,:,S2));%,Vol3(:,:,S3));
        stack_mask=cat(3,Mask1(:,:,S1),Mask2(:,:,S2));%,Mask3(:,:,S3));


        f1=[feat_Vol1(:,:,:,fselect(j))];
        f2=[feat_Vol2(:,:,:,fselect(j))];
        %  f3=[feat_Vol3(:,:,:,fselect(j))];


        stack_feat=cat(3,f1,f2);%,f3);

    end
    stack_feat=stack_feat.*stack_mask;
    stack_img=rescale_range(stack_img,0,1);
    stack_feat=rescale_range(stack_feat,0,1);


    pt=0;


    %     my4DMap(stack_img(:,:,(1:3)),stack_feat(:,:,(1:3)),stack_mask(:,:,(1:3)));
    %
    %     savefig(gcf,[a,'heatmap-',num2str(fselect(j)),sufix,'.fig'])
    %
    %
    %     my4DMap(stack_img(:,:,(4:6)),stack_feat(:,:,(4:6)),stack_mask(:,:,(4:6)));
    %
    %     savefig(gcf,[b,'heatmap-',num2str(fselect(j)),sufix,'.fig'])
    %
    %     my4DMap(stack_img(:,:,(7:9)),stack_feat(:,:,(7:9)),stack_mask(:,:,(7:9)));
    %
    %     savefig(gcf,[b,'heatmap-',num2str(fselect(j)),sufix,'.fig'])

    my4DMap(stack_img(:,:,:),stack_feat(:,:,:),cat(3,Mask1(:,:,S1),Mask2(:,:,S2)));
    vv(stack_img,stack_mask)
    savefig(gcf,[a,'-',b,'heatmap-',num2str(fselect(j)),sufix,'.fig'])


    % for i=[3,5,8]
    %
    %                 regions = stack_mask(:,:,i);
    %                 image=stack_img(:,:,i);
    %                 size(regions)
    %                 size(image)
    %
    %                 [dx,dy]=gradient(double(regions));
    %                 edgei=sqrt(dx.^2+dy.^2);
    %                 edgei=edgei~=0;
    %                 checkersegbound=rgbmaskrgb(rescale(image),edgei,[0 0 1]);
    %
    %                 size(checkersegbound)
    %                 figure
    %                 A= imshow(checkersegbound);
    %
    % end


end
end

%GUI objects
function my4DMap(vol,map_vol,mask)

figure('Color','white');

begin_index = find(mask==1,1,'first');
end_index = find(mask==1,1,'last');

[~,~,begin_slice] = ind2sub(size(mask),begin_index);
[~,~,end_slice] = ind2sub(size(mask),end_index);
slice = round((end_slice+begin_slice)/2);

slider = uicontrol('Style', 'slider',...
    'Max',size(vol,3),'Min',1,...
    'Units', 'normalized', ...
    'Position', [.25 .005 .4 .04],...
    'SliderStep', [1/(size(vol,3)-1) 1/(size(vol,3)-1)], ...
    'Value', slice, ...
    'Callback', @move_slider);

text_box = uicontrol('Style', 'text',...
    'Units', 'normalized', ...
    'Position', [.675 .006 .04 .04],...
    'String', num2str(slice));

set(gcf,'UserData',map_vol);
showmap(slice, vol, map_vol, mask);



%callback to slider
    function move_slider(~,~)
        mapdata = get(gcf,'UserData');
        slice = round(get(gcbo,'Value'));
        set(gcbo,'Value',slice);
        set(text_box,'String',num2str(slice));
        save_xlim = xlim;
        save_ylim = ylim;
        showmap(slice,vol,mapdata,mask);
        xlim(save_xlim);
        ylim(save_ylim);
    end

end

%where the images are displayed
function showmap(slice,origdata,mapdata,mask)

%----original-----%
%     subplot(1,2,1);
%      maskslice = transpose(mask(:,:,slice));
%     [dx,dy] = gradient(maskslice);
%     edgei=sqrt(dx.^2+dy.^2);
%     edgei=edgei~=0;
%     imagesc(rgbmaskrgb(rescale(transpose(origdata(:,:,slice))),edgei,[0 1 0]));colormap gray; axis off;
%     axis fill
%     %     imagesc(transpose(origdata(:,:,slice)));colormap(gray);axis off; %previously used
%     title([ name,' slice ' num2str(slice) ' of ' num2str(size(origdata,3))])

%----overlay------%
%     subplot(1,2,2);
baseslice = transpose(origdata(:,:,slice));
heatmapslice = transpose(mapdata(:,:,slice));

baseslice = baseslice/max(baseslice(:));
rgbslice = baseslice(:,:,[1 1 1]);

maskslice = transpose(mask(:,:,slice));
bwROIlocations = (maskslice ~= 0);
g = imagesc(heatmapslice);colormap(gca,'jet');
%caxis([0 1]);
colorbar;
%     g = imagesc(heatmapslice);colormap(gca,'gray'); %never again.

alpha(g,1);
hold on; h = imagesc(rgbslice);
set(h,'AlphaData',~bwROIlocations);
axis off
axis fill

axis off;
%colorbar; %optional
end

