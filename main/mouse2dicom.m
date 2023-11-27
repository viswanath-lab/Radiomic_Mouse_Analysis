%mouseDicomToMHA Convert images in an easier format to work with that
%matches the labels
clear

addpath(genpath('D:\Radiomics_Mouse_Analysis\image_functions'))

dirPath='D:\Radiomics_Mouse_Analysis\IBD_Mouse_MRIs\Dicom';
cd('D:\Radiomics_Mouse_Analysis\IBD_Mouse_MRIs\MHA')

Path{1}='D:\Radiomics_Mouse_Analysis\IBD_Mouse_MRIs\Dicom\';%Path to data one for each specific subpath
for i=1
    cd(Path{i});
    temp=dir(Path{i});
    
    for j=3:length(temp)
        mousefile=temp(j).name;
        temp2=dir([Path{i},'\',mousefile]);
        fcount=1;
        for k=3:length(temp2)
            view=temp2(k).name;
            temp3=dir([Path{i},'\',mousefile,'\',view]);
            if length(temp3)>10
                count=1;
                for p=3:length(temp3)
                    slice=temp3(p).name;
                    string=[Path{i},'\',mousefile,'\',view,'\',slice];                 
                    if (strncmp(slice,'MRIm',4))
                    string=[Path{i},'\',mousefile,'\',view,'\',slice];
                    Dinfo=dicominfo(string);
                    PD=[Dinfo.PixelSpacing;1];
                    X=double(dicomread(string));

                    Vol(:,:,count) = transpose(X);
       
                    count=count+1;
                    else
           
                    end
                end
                fname=strrep(mousefile,'DICOMS','MHA');
                fname=[fname,'_',num2str(fcount),'.mha'];
                mha_write_volume(fname,Vol,PD)
                fcount=fcount+1;
                clear Vol PD
            else
            end
            
            %             offest=info.ImagePositionPatient;
            %             fname=strrep(mousefile,'DICOMS','MHA');
        end
    end
end
