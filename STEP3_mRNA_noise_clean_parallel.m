%% STEP3 clean up the mRNA results
% the function will auto run all the data inside result folder
% the result will saved in the subfolder Noise_deleted as the same name 

close all
clear all
addpath('./function/');
addpath('./result3D/');

mat = dir('result3D/*.mat'); 
if ~exist('result3D/Noise_deleted', 'dir')
   mkdir('result3D/Noise_deleted')
end

pwd
% parpool(12)

% change q to chose specific data if needed
for q = 1:length(mat) 
    
    noise_clean(q,mat);

 end

function noise_clean(q,mat)

    disp(['clean embryo ' mat(q).name ' which is number ' num2str(q)]);
    if ~exist(['result3D/Noise_deleted/' mat(q).name], 'file')  
        load(['result3D/' mat(q).name]); 

        index = knn_filter(total_coord(:,1:3),5,20);
        in_range = total_coord((index~=0),:);
        out_range = total_coord((index==0),:);

        figure
        set(gcf, 'Position',  [100, 100, 1000, 400])
        %sgtitle(mat(q).name)
        subplot(1,3,1)
        scatter3(total_coord(:,1),total_coord(:,2),total_coord(:,3),[],total_coord(:,7),'b');
        title(mat(q).name)
        axis equal 

        subplot(1,3,2)
        title(['Inrange = ' num2str(length(in_range)) ' Outrange = ' num2str(length(out_range))])
        hold on
        scatter3(in_range(:,1),in_range(:,2),in_range(:,3),'b')
        scatter3(out_range(:,1),out_range(:,2),out_range(:,3),'r')
        axis equal

        subplot(1,3,3)
        title('cleaned noise data')
        hold on
        scatter3(in_range(:,1),in_range(:,2),in_range(:,3),'b')
        axis equal
        savefig(['result3D/Noise_deleted//' mat(q).name '.fig'])
        total_coord_cleaned = in_range;
        save(['result3D/Noise_deleted/' mat(q).name],'total_coord_cleaned')
    else
        return
    end
    
end