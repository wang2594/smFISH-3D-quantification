%% this code is used for big mRNA spots(connected individual mRNA) seperation
clear all
close all
addpath('./function/');
currentFolder = pwd;
addpath(genpath(currentFolder));
%% input 
load('mrna_runtable.mat')
load('bmp_57_18_bmp2b_03D.mat') % mRNA coordinates
load('bmp_57_18_bmp2b_lab.mat') % mRNA mask
file_n=9;
image_folder1 = ['mRNA_tiff/' char(runinglist.folder1(file_n))];
image_folder2 = char(runinglist.folder2(file_n));
foldername= char(runinglist.folder3(file_n));
file_c=char(runinglist.output_name(file_n));
last_slide=runinglist.slides_number(file_n);
z_micro=runinglist.z_stack(file_n);
%% set up parameters for the images
dropoff_scale = -0.18;
channel=3;
area_mean=mean(total_coord(:,4));
int_sum_mean=mean(total_coord(:,7));
im_total =[];
a=size(lab,1);
b=size(lab,2);
c=size(lab,3);
slice_direction = 2;
for slides = 1: last_slide
    if last_slide<100

        if slides>9
        folder_path = fullfile(image_folder1,image_folder2,foldername,[foldername '_z' num2str(slides) 'c' num2str(channel) '.tif']);
        ims=imread(folder_path);  
        else
        folder_path = fullfile(image_folder1,image_folder2,foldername,[foldername '_z0' num2str(slides) 'c' num2str(channel) '.tif']);
        ims=imread(folder_path);
        end
    else
        if slides>99
        folder_path = fullfile(image_folder1,image_folder2,foldername,[foldername '_z' num2str(slides) 'c' num2str(channel) '.tif']);
        ims=imread(folder_path);    
        elseif slides>9
        folder_path = fullfile(image_folder1,image_folder2,foldername,[foldername '_z0' num2str(slides) 'c' num2str(channel) '.tif']);
        ims=imread(folder_path); 
        else
        folder_path = fullfile(image_folder1,image_folder2,foldername,[foldername '_z00' num2str(slides) 'c' num2str(channel) '.tif']);
        ims=imread(folder_path);
        end
    end        
%%  Converts the image RGB to the grayscale 
    im=rgb2gray(ims);
%%  Convert to double.
    im = double(im);
%% scale the drop off
    if slice_direction==1
        z_position = (slides-1) * z_micro;
    elseif slice_direction==2
        z_position = (last_slide-slides+1) * z_micro;
    end   
    im(im>0) = im(im>0) - z_position*dropoff_scale;
%% Apply the filter
    H = -fspecial('log',15, 1.5);
   % Here, we amplify the signal by making the filter "3-D"
    H = 1/3*cat(3,H,H,H);
    outims = imfilter(im,H,'replicate');
    % Set all negative values to zero
    outims(outims<0) = 0;
    % Normalize ims2
    ims2 = outims/max(outims(:)); % make it from 0-1
    % stack the sections
    im_total = cat(3,im_total,ims2);
end
%% look for potential connected individual mRNA
bigspots=find(total_coord(:,4)>area_mean*1.5 & total_coord(:,7)>int_sum_mean*1.5);
n=size(bigspots,1);
lab=double(lab);
 for j=1:n
ind=find(lab==bigspots(j)); % 3D mask number and corresponding bigspots number
[row,col,ver] = ind2sub(size(lab),ind);
m=size(row,1);
onemRNAspot=zeros(a,b,c);
for i=1:m
    onemRNAspot(row(i),col(i),ver(i))=im_total(row(i),col(i),ver(i));
end
% extract single mRNA spot and find the local maximum in the surrounding area
% area_mean=mean(total_coord(:,4))
% r=(area_mean*3/(4*pi)).^(1/3)
[subs,vals] = nonMaxSupr( onemRNAspot, [2 2 2], 0.1, [] );% set the radius based on the averaged area
total_coord(bigspots(j),10)=size(subs,1); % number of local maximum
Pos_max{bigspots(j),1}=subs; % save the location 
Pos_max{bigspots(j),2}=vals; % save the maximum intensity
 end
 save('result3D/bmp_57_18_3Dsep.mat','total_coord','Pos_max', '-v7.3');