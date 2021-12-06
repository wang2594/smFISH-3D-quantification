%% this code is used for finding the nascent mRNA
% we consider the mRNA spots is nacent mRNA if it locates in the nuclear
% and sum of the intensity, area are higher than 1.5 times of the averaged
% level. Also, maximum number of nascent mRNA in each nuclear is 4.
clear
load('bmp_57_18_3Dsep.mat') % mRNA coordinates
load('nuclearmask_18_above1.mat') % nuclear masks: run SEG3D for nuclear segmentation,>>load('stack2_t1.mat')>> mask_15_nuclear=seg_result.cf_mark3d;
format bank
total_coord=total_coord_cleaned;
% x, y, z normalized to pixels to compared with nuclear mask
total_coord(:,1)=round(total_coord(:,1)./0.312); 
total_coord(:,2)=round(total_coord(:,2)./0.312);
total_coord(:,3)=round(total_coord(:,3)./3);
zmax=max(total_coord(:,3));
sum_avein=mean(total_coord(:,7));
mean_area=mean(total_coord(:,4));
for i= 1:zmax
    mask=mask_18_nuclear(:,:,i); % need to be changed mannually
    a=find(total_coord(:,3)==i);
    mRNAspots=total_coord(a,:);
    num_spot=size(mRNAspots,1);
    for j=1:num_spot
    if (mask(mRNAspots(j,1),mRNAspots(j,2))>0 && total_coord(a(j),7)>1.5*sum_avein && total_coord(a(j),4)>1.5*mean_area )
       total_coord(a(j),8)=1; % 1 represent this mRNA is nascent mRNA
       total_coord(a(j),9)=mask(mRNAspots(j,1),mRNAspots(j,2)); %get the nuclear mask number of nacent mRNA location at each z section 
    end
    end
end
%% maximum number of nascent mRNA in each nuclear is 4
s=size(total_coord,1);
nasent_number=zeros(s,1);
for i = 1:s
    if total_coord(i,9)>0
    spots_rem=find(total_coord(:,9) ==total_coord(i,9));
    n=size(spots_rem,1); % number of mRNA spots in the nuclear
    nasent_number(i,1)=n;
    if n>=5
    [sub_sort,index]=sortrows(total_coord(spots_rem,:),7, 'descend');% mRNA spots ranking by the sum of the intensity
    re=spots_rem(index(5:n),1);
    total_coord(re,8)=0;
    total_coord(re,9)=0;
    end
    end
end 
%% mature and nascent mRNA intensity
total_coord(total_coord(:,10)==0,10)=1; %single mRNA spots with smaller size and intensity are also single individual mRNAs
total_coord(:,11)=total_coord(:,7)./total_coord(:,10); % get the intensity of individual mRNA
Ave_int_maturemRNA=mean(total_coord(total_coord(:,8)==0,11)); % averaged intensity of mature mRNA
total_coord(total_coord(:,8)==1,10)=round(total_coord(total_coord(:,8)==1,7)/(Ave_int_maturemRNA*0.5));% quantification of nascent mRNA
total_coord(total_coord(:,8)==1,11)=total_coord(total_coord(:,8)==1,7)./total_coord(total_coord(:,8)==1,10);% quantification of nascent mRNA
%% total_coord = [1_x,2_y,3_z,4_area pixal number,5_max_intansity, 6_avg_intansity,7_sum_intansity,8_nascent mRNA(1), 9_nuclear mask number,10_individual mRNA number of mature and nascent mRNA,11_intensity of individual mRNA]
%% x,y,z scaled
total_coord(:,1)=total_coord(:,1)*0.312;% x: 0.312 um for each pixel
total_coord(:,2)=total_coord(:,2)*0.312;% y: 0.312 um for each pixel
total_coord(:,3)=total_coord(:,3)*3;% z: 3 um for each pixel
total_coord_cleaned=total_coord;
%save('bmp_57_18_bmp2b_cleaned','total_coord_cleaned') %change saving name