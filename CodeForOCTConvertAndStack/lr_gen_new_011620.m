clc, clear, close all
tic

for index = 23:23
    name = ['Default_00' num2str(floor(index/10)) num2str(mod(index,10)) '_Mode3D.oct'];
    path = [   ''];
    file_name = [path name];
    handle = OCTFileOpen(file_name);

    %%%%%%%%%%%%%%
     disp(['-----------------index:' num2str(index) '--------------']);
     toc
    %%%%% reading chirp vector %%%%%%

    Chirp = OCTFileGetChirp(handle);
 
    %%%%% reading spectral raw data %%%%%%
%%
    
    NrRawData = OCTFileGetNrRawData(handle);
%     spectrumData3D = zeros(2048,1500,NrRawData);
%     windata3D = zeros(2048,1500,NrRawData);
    %%
    
    
    for ind = 1:10:NrRawData
        
        if mod(ind,25) == 0          
            disp(['------------finish: '  num2str(ind/NrRawData) '--------------']);
            toc
        end
        
        [RawData, Spectrum] = OCTFileGetRawData(handle, ind-1);    
        Background  =  repmat(Spectrum, [1 size(RawData,2)]);
        [Chirp ] = OCTFileGetChirp( handle );    
        SubtractedData=RawData-Spectrum;%already subtracted in OCTFileGetRawData.    
        spectrumData3D(:,:,ind) = SubtractedData;
             
        
        %fully reconstructed
        %SubtractedData=RawData;
        Apod = SubtractedData;
        dims = size(Apod);
        
        %spectrum with 1/4 spectrum
        start_hav = round(dims(1)/2) - round(dims(1)/8);
        end_hav = start_hav + round(dims(1)/4) - 1;
        
        Apod2 = zeros(dims(1), dims(2));
        Apod2(start_hav:end_hav, :) = Apod(start_hav:end_hav, :,1);
        
        %spectrum with 1/8 
        start_hav = round(dims(1)/2) - round(dims(1)/16);
        end_hav = start_hav + round(dims(1)/8) - 1;
        
        Apod3 = zeros(dims(1), dims(2));
        Apod3(start_hav:end_hav, :) = Apod(start_hav:end_hav, :,1);
        
        if numel(dims) == 2
              for j = 1:dims(2)
                k(:,j) = lamb2k_v2(squeeze(Apod(:,j)), Chirp); %individual Aline
                k2(:,j) = lamb2k_v2(squeeze(Apod2(:,j)), Chirp); 
                k3(:,j) = lamb2k_v2(squeeze(Apod3(:,j)), Chirp);       
              end
             data = k;
             data2 = k2;
             data3 = k3;

        else
            h = waitbar(0,'Preparing raw volume data for processing...');
            %3 Dimensional Data
            for n = 1:dims(3)
                for j = 1:dims(2)
                    k(:,j,n) = lamb2k_v2(squeeze(Apod(:,j,n)), Chirp); %individual Aline
                 %   T(:,j,n)=fft(lamb2k_v2(squeeze(Apod(:,j,n))), Chirp);
                end
                waitbar(n/dims(3),h,[num2str(dims(3)-n) ' frames remaining for interpolation']);
            end
            data = k;
        end

        [x,y]=size(data);
        Windata=zeros(x,y);
        w = hann(x) ;
        for i=1:y
            Windata(:,i) = w.*data(:,i);
        end
        data = Windata(:,1:end);
%         windata3D(:,:,ind)=Windata;
        A = fft(data);
        A = A(1:1024,:,:);
        Spectraldata = data;
        ComplexData = A;         

        imOut = 10*log10(abs(ComplexData));
%         RawDataNFrame(:,:,ind) = data;
%        imOut = imadjust(uint8(imageOut), [30 55]/255);
        
        windowWidth = 3;
        kernel = ones(windowWidth)/(windowWidth*windowWidth);
%         imOut = uint8(conv2(uint8(imOut), kernel, 'same')/(windowWidth*windowWidth));
        
        figplot = 0;
        if figplot == 1
         figure(1), imagesc(imOut);colormap(gray)
        end
        %
        w2 = hann(x) ;
        for i=1:y
            Windata2(:,i) = w2.*data2(:,i);
        end
        data2 = Windata2(:,1:end);
%         windata3D(:,:,ind)=Windata;
        A2 = fft(data2);
        A2 = A2(1:1024,:,:);
        Spectraldata = data2;
        ComplexData2 = A2;         
%         RawDataNFrame(:,:,ind) = data;
        imOut2 = 10*log10(abs(ComplexData2));
%         imOut2 = uint8(conv2(uint8(imOut2), kernel, 'same')/(windowWidth*windowWidth));
        
        
        if figplot == 1
            figure(2), imagesc(imOut2);colormap(gray)
        end
        
        
        %
        w2 = hann(x) ;
        for i=1:y
            Windata2(:,i) = w2.*data3(:,i);
        end
        data3 = Windata2(:,1:end);
%         windata3D(:,:,ind)=Windata;
        A3 = fft(data3);
        A3 = A3(1:1024,:,:);
        Spectraldata = data3;
        ComplexData3 = A3;         
%         RawDataNFrame(:,:,ind) = data;
        imOut3 = 10*log10(abs(ComplexData3));
%         imOut3 = imadjust(uint8(imageOut3), [30 55]/255);
        
        if figplot == 1
            figure(3), imagesc(imOut3);colormap(gray)
        end
        HRx1(:,:,ind) = imOut(1:512, :);
        LRx4(:,:,ind) = imOut2(1:512, :);
        LRx8(:,:,ind) = imOut3(1:512, :);
    end  
    
    
%%
%     figure(4);clf;
%     plot(Spectrum);
%     figure(5);clf;
%     imagesc(RawData);

    
    

%     data = spectrumData3D; 
% %     save SpectrumData data
%     out_name = [path '\Nframe_Complex_Mode_' num2str(index)  '.mat'];  
% % 
%     save(out_name, 'windata3D' ,'spectrumData3D', 'Background','-v7.3')
% %     save('ComplexData2')

    out_name = [path '' num2str(index)  '.mat'];  
% 
    save(out_name, 'HRx1' ,'LRx4', 'LRx8','-v7.3')
    
    toc
    clear Hrx1
    clear LRx4
    clear LRx8 
    
    
end