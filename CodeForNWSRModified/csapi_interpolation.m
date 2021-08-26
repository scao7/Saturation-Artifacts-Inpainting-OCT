clc;
clear all;
close all;
f=imread('./outs/W5NWSR_1xInterp/HH600averaged.tif'); 
ff=imread('./outs/W5NWSR_1xInterp/out_onion_LL276_CSI.tif'); 
[R,C]=size(ff);
emtyVector = zeros(R,1);
emtyVector(emtyVector ==0) = 255;
experiment_name='NWSR_1xInterp';
outfolder=['outs/W5' experiment_name];
outfname=[outfolder,'/out_onion_LL276_CSI_SR'];
missingcol=[];
findcnt=1;
for ii=1:C
    colidx=zeros(1,4);
    value=zeros(R,4);
    
    if ff(:,ii) == emtyVector
        missingcol(1,findcnt)=ii;
        fprintf('%d: %d\n',findcnt,ii);
        findcnt=findcnt+1;
        idx=1;
        count=1;
        while count<=2
            if ff(:,ii-idx) ~=emtyVector
                colidx(1,idx)=ii-idx;
                value(:,idx)=ff(:,ii-idx);
                idx=idx+1;
                count=count+1;
            else 
                idx=idx+1;
            end
        end
        idx=1;
        count=1;
        while count<=2
            if ff(:,ii+idx) ~=emtyVector
                colidx(1,2+count)=ii+idx;
                value(:,2+count)=ff(:,ii+idx);
                idx=idx+1;
                count=count+1;
            else 
                idx=idx+1;
            end
        end
    
        columnvalue=cubicinterp(colidx,value,ii,R);
        ff(:,ii)=columnvalue;
    end


end

figure(1)
imagesc(ff)
colormap('gray')
imwrite(ff,[outfname '.tif'],'tif');
figure(2)
imagesc(f)
colormap('gray')
 %     fnplt(cs);
%     hold on
%     plot(X,Y,'o')
%     legend('cubic spline','data')
%     hold off