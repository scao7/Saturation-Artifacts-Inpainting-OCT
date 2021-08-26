%双三次插值具体实现
clc,clear;
ff=imread('./datasets/ourDatasets/CreatedImages/without23/HH600averaged.tif'); 
% ff =rgb2gray(fff);%转化为灰度图像
% ff(:,1)=inputimage(:,749);
% ff(:,2)=inputimage(:,751);
[mm,nn]=size(ff);               %将图像隔行隔列抽取元素，得到缩小的图像f
m=mm/2;
n=nn/2;
% f=ff;
f =zeros(m,n);
for i=1:m
   for j=1:n
     f(i,j)=ff(2*i,2*j);
   end
end
 
k=5;                       %设置放大倍数
bijiao1 =imresize(f,k,'bilinear');%双线性插值结果比较
bijiao =uint8(bijiao1);
 
a=f(1,:);
c=f(m,:);             %将待插值图像矩阵前后各扩展两行两列,共扩展四行四列
b=[f(1,1),f(1,1),f(:,1)',f(m,1),f(m,1)];
d=[f(1,n),f(1,n),f(:,n)',f(m,n),f(m,n)];
a1=[a;a;f;c;c];
b1=[b;b;a1';d;d];
ffff=b1';
f1=double(ffff);
g1 =zeros(k*m,k*n);
for i=1:k*m                 %利用双三次插值公式对新图象所有像素赋值
    u=rem(i,k)/k;
    i1=floor(i/k)+2;
    A=[sw(1+u) sw(u) sw(1-u) sw(2-u)];  
    for j=1:k*n
        v=rem(j,k)/k;
        j1=floor(j/k)+2;
        C=[sw(1+v);sw(v);sw(1-v);sw(2-v)];
        B=[f1(i1-1,j1-1) f1(i1-1,j1) f1(i1-1,j1+1) f1(i1-1,j1+2)
            f1(i1,j1-1)   f1(i1,j1)  f1(i1,j1+1)   f1(i1,j1+2)
            f1(i1+1,j1-1)   f1(i1+1,j1) f1(i1+1,j1+1) f1(i1+1,j1+2)
            f1(i1+2,j1-1) f1(i1+2,j1) f1(i1+2,j1+1) f1(i1+2,j1+2)];
        g1(i,j)=(A*B*C);
    end
end
g=uint8(g1); 
 
imshow(uint8(f));title('缩小的图像');             %显示缩小的图像
figure,imshow(ff);title('原图');               %显示原图像
figure,imshow(g);title('双三次插值放大的图像');     %显示插值后的图像
figure,imshow(bijiao);title('双线性插值放大结果');     %显示插值后的图像 
% mse=0;
% ff=double(ff);
% g=double(g);            
% ff2=fftshift(fft2(ff));   %计算原图像和插值图像的傅立叶幅度谱                            
% g2=fftshift(fft2(g));
% figure,subplot(1,2,1),imshow(log(abs(ff2)),[8,10]);title('原图像的傅立叶幅度谱');
% subplot(1,2,2),imshow(log(abs(g2)),[8,10]);title('双三次插值图像的傅立叶幅度谱');
im =ff;
im_out =g;
PSNR=comp_psnr(im,im_out);
fprintf(' PSNR = %d\n',PSNR); 
% 基函数代码：
function [PSNR,SSIM,ISNR]=comp_psnr(im,imf,imn)
% MSE=mean(mean((im-imf).^2));
    [rows columns] = size(im);
    squaredErrorImage = (double(imf) - double(im)) .^ 2;
    MSE= sum(squaredErrorImage(:)) / (rows * columns);

    if max(im(:))<2
        MaxI=1;
    else
        MaxI=255;
    end
    PSNR=10*log10((MaxI^2)/MSE);
    %
    if nargout>1
        SSIM=ssim_index(im,imf);
    end
    if nargout>2 && (exist('imn','var') || ~isempty(imn))
        ISNR=10*log10((mean(mean((im-imn).^2))/MSE));
    end
end

function A=sw(w1)
    w=abs(w1);
    if w<1&&w>=0
        A=1-2*w^2+w^3;
    elseif  w>=1&&w<2
        A=4-8*w+5*w^2-w^3;
    else
        A=0;
    end
end