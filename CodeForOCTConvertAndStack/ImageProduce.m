clc, clear, close all
tic

for index = 10:10  
    name = ['Default_00' num2str(floor(index/10)) num2str(mod(index,10)) '_Mode3D.oct'];
    path = ['../'];
    file_name = [path name];
    handle = OCTFileOpen(file_name);
    disp(['-----------------index:' num2str(index) '--------------']);
    toc
    Chirp = OCTFileGetChirp(handle);
    NrRawData = OCTFileGetNrRawData(handle);
    Intensity = OCTFileGetIntensity(handle);
    convertedImage = uint8(Intensity);
   
    
    for order = 1:NrRawData
        
        ind = order;
%         if order == 1
%             ind = 222;
%         end
%        
%         if order ==2
%             ind = 403;
%         end
%         
%         if order == 3
%             ind = 648
%         end
%         
        if mod(ind,25) == 0          
            disp(['------------finish: '  num2str(ind/NrRawData) '--------------']);
            toc
        end

        [RawData, Spectrum] = OCTFileGetRawData(handle, ind-1);     

        
        Background  =  repmat(Spectrum, [1 size(RawData,2)]);
        [Chirp ] = OCTFileGetChirp( handle );    
        SubtractedData=RawData;%already subtracted in OCTFileGetRawData.    
        spectrumData3D(:,:,order) = SubtractedData;
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
        im_temp = imadjust(convertedImage(:,:,order), [5 80]/255);
        imOutIntensity(:,:,order) = im_temp;
        
        imOut = 10*log10(abs(ComplexData));
        imOut3D(:,:,order) = imOut;
    end
end
