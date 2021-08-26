%%
testImg_indices=1:20;
scale_factor=1;
experiment_name='NWSR_1xInterp'
%
NN=length(testImg_indices);% number of test images
%
% pre-allocating memory for creating a table of measures.
qm_data=zeros(NN,2);% #images by #quality measures & time
row_name=cell(NN,1);
%
addpath('C:\Program Files\MATLAB\R2020b\toolbox\ompbox10');

for i=1
    strnumber=num2str(testImg_indices(i));
    fprintf('image #%s\n',strnumber)
    pth=['./outs/W5NWSR_1xInterp']; 
    % ** Load a test and its HH image
    testfile=['out_onion_LL276_CSI.tif'];
%     testfile=['HH600experiment',strnumber,'.tif'];
    cleanfile='LL276.tif'; 
    imn = single(imread(fullfile(pth,testfile)));
    im= single(imread(fullfile(pth,cleanfile)));
    % the outputs are saved as images in the following folder
    outfolder=['outs/W5' experiment_name];
    if ~exist(outfolder,'dir')
        mkdir(outfolder);
    end
    NUM=i;
    NUM_neighbors=20;
    outfname=[outfolder,'/out_onion_LL276_CSI_SR',strnumber];% output file name
    % 
    ps=[8 8]; % patch size for sparse representation
    psnlm=[6 6];% patch size for NLM filtering
    dlfile='dictionary_onion'; % dictionary file name
   
    % run the algorithm
    [im_out,time_end]=main_reconstruct_oct_image(i,imn,...
        scale_factor,ps,psnlm,dlfile,NUM,NUM_neighbors);
    
    PSNR=comp_psnr(im,im_out);
    fprintf('i%d: %d\n',i, PSNR);
    qm_data(i,:)=[PSNR,time_end];
    row_name{i}=strnumber;
    %
    imwrite(im2uint8(im_out/255),[outfname '.tif'],'tif');
end
% % Draw a table
% f=figure;
% colname={'PSNR','TIME'};
% colformat=repmat({'short g'},1,numel(colname));
% t = uitable('Parent', f,'Data', qm_data,'RowName',row_name,...
%     'ColumnName', colname,'ColumnFormat', colformat,...
%     'Units','normalized','Position',[0 0 1 1]);
% t.FontSize=12;
% f.Name=[ experiment_name ' - 18 foveal images'];
