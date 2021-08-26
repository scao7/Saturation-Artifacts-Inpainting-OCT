% clc,clear,close all
% 
% load('imOut3D.mat' );
% load('spectrum3D.mat');
% load('Chirp.mat')


for imageIndex = 1:750
    BscanOriginal  = imOutIntensity(:,:,imageIndex);
    maxAdjust = max(max(BscanOriginal));
    minAdjust = min(min(BscanOriginal));
%     BscanOriginal = imadjust(BscanOriginal, [5 80]/255);
% 
    f1 = figure(1) 
    imshow(BscanOriginal)

    %% find the artifacts 
    Spectrum = spectrumData3D(:,:,imageIndex);
    peak =[];
    for x = 1:1500
        peak(1,x) = max(Spectrum(:,x));   
    end

    index =1 ;
    columnWithArtifacts  = [];
    for x = 1:1500
       if peak(1,x) >=  9.996508993812499e+04
          columnWithArtifacts(index) = x
          index = index +1;
       end
    end

    BscanMask = BscanOriginal;
    BscanMask(:,columnWithArtifacts)=255;
    f2 = figure(2)
    imshow(BscanMask) 
    combImg=imfuse(BscanOriginal,BscanMask,'montage');
%     f3 = figure(3)
%     imshow(combImg);
    combImg3D(:,:,imageIndex) = combImg;  
end


