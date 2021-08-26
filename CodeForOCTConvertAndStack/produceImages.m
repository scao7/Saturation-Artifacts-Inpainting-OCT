clc,close all,clear;
load('23.mat');
% figure(1)
% imshow(uint8(HRx1(:,:,1))*2);
% figure(2)
% imshow(uint8(LRx4(:,:,1))*2);

for i= 1:10:740
    
   imwrite(uint8(HRx1(:,:,i))*2, [ 'HRx1_optical/' '23_' num2str(i) '.png']);
   imwrite(uint8(LRx4(:,:,i))*2, [ 'LRx4_optical/' '23_' num2str(i) '.png']);
   
end