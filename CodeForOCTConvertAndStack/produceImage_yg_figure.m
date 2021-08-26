close all, clc
% load('Default_0023_Mode3D.mat');
count = 1;
train_set = 4:2:23;
test_set = 3:2:23;
for ii = 10:10
tic 
ii
index = ii;
name = ['HR_LR_pairs_' num2str(index) '.mat'];
load(name);
    win = 3;
%%
f = filesep;
for i = 450:450
    image_hr = uint8(mean(HRx1(:,:,i:i+4),3));
    image_hr = uint8(conv2(double(image_hr),ones(win)/(win*win), 'same'));
%     figure(1), imagesc(image_hr), colormap(gray)
    image_hr = imadjust(uint8(image_hr), [120 230]/255);
  
%     figure(2), imagesc(image_hr2), colormap(gray)
    
    image_hr = cat(3,image_hr,image_hr,image_hr);
    imwrite(image_hr,'Full_hr.png');


    %x4

    image_lr_4 = uint8(mean(LRx4(:,:,i:i+4),3));
%     figure(1), imagesc(image_lr_4), colormap(gray);
    image_lr_4 = conv2(double(image_lr_4),ones(win)/(win*win), 'same');

    image_lr_4 = imadjust(uint8(image_lr_4), [120 230]/255);

%      figure(2), imagesc(image_lr_4), colormap(gray);
%     image_lr_4 = uint8(image_lr_4(1:4:end,1:4:end));
%      figure(3), imagesc(image_lr_4), colormap(gray);
    image_lr_4 = cat(3,image_lr_4,image_lr_4,image_lr_4);
    imwrite(image_lr_4,'full_lr.png');
   
    lr_down_4 = image_lr_4(1:4:end, 1:4:end,:);
    lr_down_4 = replem(lr_down_4, 4,4);
    imwrite(lr_down_4,'full_lr_downsample.png');
    
    %x8
%     image_lr_8 = LRx8(1:8:end,1:8:end,i);
%     image_lr_8 = cat(3,image_lr_8,image_lr_8,image_lr_8);
%     imwrite(image_lr_8,['imageNew'  f 'test_v2_lrx8' f  int2str(count) '_lrx8.png']);

    count = count + 1;
    toc
end

clear spectrumData3D;

end