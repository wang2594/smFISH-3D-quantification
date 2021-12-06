%% this code is used for mRNA spots 3D segmentation from the tif files of 
% whole mount embryos.
% example image "Wholemount_RNA_57hpf_PBST_18-Stitching_z01c1": 
%   01: z-stack section; 
%   c1: channel 1(Chd); 
%   c2: channel 2(DAPI);
%   c3: channel 3(bmp2b)
% All the input information are saved in the 'mrna_runtable.mat' file.
% you can add your own information into the table.
% Images of whole mount embryos are segmentated by mRNA_quant_3D

clear all
close all
addpath('./function/');
pwd
% parpool(12); % parallel running if number of file_n is more than 1

load('mrna_runtable.mat')

for file_n =9
    file_n
    run_mRNA(file_n,runinglist);
end

function run_mRNA(file_n,runinglist)
%% set up parameters for the images
    threshold =0.071;  % threshold for smallest intansity
    dropoff_scale = -0.18; % scale the drop off  im = im- z_position*dropoff_scale;
    slice_direction = 2; % set up the slice (animal pole -vege pole)  == 1  // (vege pole- animal pole) == 2
                         % check if the first slice is on the top of the 3D
    use_filter =1;  % check if use gauss filter or not 1== use // 2== not use
%% input 
    image_folder1 = ['mRNA_tiff/' char(runinglist.folder1(file_n))];
    image_folder2 = char(runinglist.folder2(file_n));
    foldername= char(runinglist.folder3(file_n));
    file_c=char(runinglist.output_name(file_n));
    last_slide=runinglist.slides_number(file_n);
    x_pixel=runinglist.xy_resolution(file_n);
    y_pixel=runinglist.xy_resolution(file_n);
    z_micro=runinglist.z_stack(file_n);
    channel1 = char(runinglist.channel1(file_n));     
    channel3 = char(runinglist.channel3(file_n));     

%% set up running channel
    for channel=3 
        
    if channel==1
        filename=['result3D/' file_c '_' channel1 '_03D.mat'];
        filename_lab=['result3D/' file_c '_' channel1 '_lab.mat'];
        split_large_spot = 'n';  
        channel_name = channel1;
    elseif channel==3
        filename=['result3D/' file_c '_' channel3 '_03D.mat'];
        filename_lab=['result3D/' file_c '_' channel3 '_lab.mat'];
        split_large_spot = 'n';
        channel_name = channel3;
    end
  
%% start processing
    %initialize coordinate matrix
    total_coord=[];
    coord = [];
%% loop for each section
    % get the single slide coordinate
    %coord=[x,y,spot size,max_intansity,avg_intansity,sum_intansity)]; x,y not sclale
    [coord,lab]=mRNA_quant_3D(image_folder1,image_folder2,foldername,channel,last_slide,threshold,dropoff_scale,z_micro,slice_direction,use_filter);
   
    % scaled to real size on micro meter
    if isempty(coord)~=1
        one_coord=[coord(:,1)*x_pixel,coord(:,2)*y_pixel,coord(:,3)*z_micro,coord(:,4),coord(:,5),coord(:,6),coord(:,7)];
        total_coord=[total_coord;one_coord];
    end
    % add to total coordinate
    % total_coord = [x,y,z,area pixal number,max_intansity, avg_intansity,sum_intansity]

    %% plot out result
    figure
    scatter3(total_coord(:,1),total_coord(:,2),total_coord(:,3),[],total_coord(:,7),'filled');
    axis equal
    colorbar
    title('total_intansity')
    %savefig(['result3D/figure/' file_c '_' channel_name '_03D.fig'])
    save(filename,'total_coord'); % save total coordinates
    save(filename_lab,'lab', '-v7.3');% save all the 3D mRNA spots segmentation mask
    end
end


