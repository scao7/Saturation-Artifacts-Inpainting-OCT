% illustrate the process from .OCT to .PNG

clc, clear, close all
tic

%%
index =23; 
name = ['Default_00' num2str(floor(index/10)) num2str(mod(index,10)) '_Mode3D.oct'];
path = [''];
file_name = [path name];

% handle = OCTFILEOPEN( filename ) Open .oct file located at filename
handle = OCTFileOpen(file_name);

disp(['-----------------index:' num2str(index) '--------------']);

toc
% data = OCTFILEGETCHIRP( handle ) Get the chirp vector in an OCT file
Chirp = OCTFileGetChirp(handle);

% OCTFILEGETNRRAWDATA  Get number of spectral raw data files in an .oct file.
NrRawData = OCTFileGetNrRawData(handle);

% OCTFILEGETINTENSITY  Get the intensity data from an .oct file
Intensity = OCTFileGetIntensity(handle);
convertedImages = uint8(Intensity);
figure(1)

% display one of the image
imshow(convertedImages(:,:,1));

% RawData, Spectrum = OCTFILEGETRAWDATA( handle, spectrumIndex, cameraIndex ) Get spectral raw data from .oct file
[RawData, Spectrum] = OCTFileGetRawData(handle, 0);
figure(2)
subplot(2,1,1)
plot(RawData)
subplot(2,1,2)
plot(Spectrum)

%% 
% this is the basic why to get specturm, depends on different images 
% some has artifacts and some don't 
% here is an exmaple for satrate artifacts image
load('spectrum3D.mat')
figure(3)
plot(spectrumData3D(:,1032,3));
% the purpose of this step is to find a threshold. 

%% the following steps is to detect artifacts
load('imOut3D.mat' );
load('spectrum3D.mat'); % this specturm3D is different from privous. This should have 750 spectrums. 
load('Chirp.mat')
for imageIndex = 1:3
    BscanOriginal  = imOut3D(:,:,imageIndex);
    maxAdjust = max(max(BscanOriginal));
    minAdjust = min(min(BscanOriginal));
    f1 = figure(4) 
    imagesc(BscanOriginal)
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
    figure(5)
    imagesc(BscanMask) 
    combImg=imfuse(BscanOriginal,BscanMask,'montage');
    figure(6)
    imshow(combImg);
    combImg3D(:,:,imageIndex) = combImg;  
end