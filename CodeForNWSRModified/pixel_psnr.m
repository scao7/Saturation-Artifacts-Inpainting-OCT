im =imread('HH600averaged.tif');
im_out =imread('out20_cubicinterpolation.tif');

%%extract im pixels corresponding to missing value
for ii=1:C
    if im(:,ii) == emtyVector
        missingcol(1,findcnt)=ii;
        fprintf('%d: %d\n',findcnt,ii);
        columnvalue(:,findcnt)=im(:,ii);
        findcnt=findcnt+1;
    end
end
im_pixel=columnvalue;

%%extract im_out pixels corresponding to missing value
for ii=1:findcnt-1
    columnvalue_out(:,ii)=im_out(:,missingcol(1,ii));
end
im_out_pixel=columnvalue_out;

PSNR=comp_psnr(im_pixel,im_out_pixel);
fprintf(' PSNR = %d\n',PSNR);

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
