

%% About
% extract 2D features from mouse data per slice
% framework can be easily adapted to other datasets


%% Setup
clear
close all;
%change to main mouse directory
currdir=mfilename('fullpath');

coreDir=currdir(1:end-49);
dataDir = [coreDir,'\IBD_Mouse_MRIs\'];%
addpath(genpath(coreDir))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAMPFILES
patients = dir([dataDir,'\SAMP_MHA\HighResAxial\']);
patients = {patients.name};
patients = patients(3:end);
volnames = patients((cellfun(@isempty,strfind(patients,'-label'))));

VOL{1}=volnames;
clear patients volnames;

patients = dir([dataDir,'\AKR_Healthy_SAMP_MHA\HighResAxial\']);
patients = {patients.name};
patients = patients(3:end-1);
volnames = patients((cellfun(@isempty,strfind(patients,'-label'))));

VOL{2}=volnames;
clear patients volnames;

patients = dir(dataDir,'\SAMP_MHA\Test\HighResAxial\');
patients = {patients.name};
patients = patients(3:end);
volnames = patients((cellfun(@isempty,strfind(patients,'-label'))));

VOL{3}=volnames;
clear patients volnames;

%Group3
addpath(dataDir,'\Post-Group4')
patients = dir([dataDir,'\Post-Group4']);
patients = {patients.name};
patients = patients(3:end);
volnames = patients((cellfun(@isempty,strfind(patients,'-label'))));
VOL{4}=volnames;
clear patients volnames;
%% Run Feature extraction Code
%for ij=1:length(VOL)
fail={}
for ij=1:1
    volnames=VOL{ij};
    clear  featstats features
    %     volnames=shit;
    for ii =1:length(volnames)
        ptvolpath = [volnames{ii}]
        ptmaskpath = [volnames{ii}(1:end-4),'-label.mha'];
        fprintf('\n Patient %s:',volnames{ii});
        if exist(ptmaskpath)==2
            %load patient data
            fprintf('\tloading data\n');
            mha_info=mha_read_header(ptvolpath)
            imgvol = double(mha_read_volume(mha_info));
            mha_info=mha_read_header(ptmaskpath)
            maskvol = double(mha_read_volume(mha_info));
            
            %find annotation slices
            [~,~,slices] = ind2sub(size(maskvol),find(maskvol==2));
            annotation_id = 2
            
            if length(slices)==0
                [~,~,slices] = ind2sub(size(maskvol),find(maskvol==4));
                annotation_id = 4
            end
            if length(slices)==0
                [~,~,slices] = ind2sub(size(maskvol),find(maskvol==1));
                annotation_id = 1
            end
            if length(slices)>0
                slices = unique(slices);
                size(slices)
                PD=mha_info.PixelDimensions();
                
                for jj = 1:length(slices)%Extract features for each slice independently
                    
                    %make 2D
                    jj
                    size(slices)
                    img = imgvol(:,:,slices(jj));
                    mask = maskvol(:,:,slices(jj));
                    %extract features
                    fprintf('\tCalling Feature Extraction Code for slice %d of %d\n',jj,length(slices));
                    
                    ws_options=[3,5,7]
                    class_options = {'gray','gradient','haralick','gabor','collage'};
                    [featints_temp,featnames,featstats_temp,statnames,feature_stack_temp] =extract2DFeatureInfo(img,mask==annotation_id,class_options,ws_options);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%for use with old samp data%%%%
                    featints_temp=featints_temp(1,(2:102));
                    featnames=featnames(1,(2:102));
                    featstats_temp=featstats_temp(:,(2:102));
                    feature_stack_temp=feature_stack_temp(:,:,(2:102));
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    bowel_wall=img.*(mask==annotation_id);
                    
                    B1=min(min(bowel_wall));
                    B2=max(max(bowel_wall));
                    
                    
                    mcount=mask(mask==annotation_id);
                    mcount=PD(1)*PD(2)*length(mcount);
                    
                    f40s={'area','area','area','area'}';
                    ffill=mcount.*ones(4,1);
                    
                    features(:,:,jj,(1:101))=feature_stack_temp;
                    featstats_temp=[featstats_temp,ffill];
                    featstats(:,:,jj) =featstats_temp;
                    featints{jj}=featints_temp;
                    features(:,:,jj,size(feature_stack_temp,3)+1)=(mask==annotation_id);
                    features(:,:,jj,size(feature_stack_temp,3)+2)=bowel_wall;
                    
                    statnames=[statnames,f40s];
                    fprintf('\tSaving featstats...');
                    %             save([volnames{ii}(1:end-4),post,'.mat'],'featstats','statnames','features','slices','imgvol','maskvol');
                end
                save([volnames{ii}(1:end-4),post,'.mat'],'featstats','statnames','features','featints','slices','imgvol','maskvol');
                clear featints_temp featnames featstats_temp statnames feature_stack_temp features featstats statnames slices imgvol maskvol features
            else
                fail{end+1}=ptvolpath;
            end
            fprintf('Done.\n');
        end
    end
    
end
clear t1 t2 f1 f2 img mask imgvol maskvol normal feat*
fprintf('***COMPLETE***\n')


%%



