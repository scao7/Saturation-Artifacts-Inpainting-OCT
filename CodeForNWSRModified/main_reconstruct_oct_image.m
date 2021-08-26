%
% imn: input image
% scale_factor: downsampling factor
% ps: patch size
% psnlm: patch size used for non-local mean denoising method
% (patchfiltering function by [1])
% dlfile: dictionary file name or path
% NUM:
% NUM_neighbors:
% im: 
% down_sample: 
%
% NOTE: if 'im' is available, the algorithm does not use off-the-shelve 
% 
% ashkan
function [im_out,time_end]=main_reconstruct_oct_image(i,imn,...
    scale_factor,ps,psnlm,dlfile,NUM,NUM_neighbors,im,down_sample)
% num = 4 or 6
% num_neighbors=15
%% set parameters
fprintf('set parameters ....\n');
if ~exist('im','var')|| isempty(im)
    im=[];
end
if ~exist('NUM','var')|| isempty(NUM)
    NUM=6;
    NUM_neighbors=15;
end
if ~exist('down_sample','var')|| isempty(down_sample)
    down_sample=1;
end
% dlfile='dicts_comp_noDC_e9_it10_4x8_odct.mat';
load(dlfile);
% if you want a dictionary with fewer columns
% [U,S,V]=svd(D{1});
% D{1}=U(:,1:63);
%
fprintf('dictionary file name: %s\n',dlfile); % for loading or saving

% number of clusters. Note that we have only one cluster (k=1).
k=1;
% ps=[4 8];
step=1;
% ** sparse coding for inpainting
sparse_func=@sparse_inp_momp;% mexOMP (denoise and inpaint) - minibatch
par_test=cell(1,k);
par_test{k}.Tdata=2; % sparsity level for cluster k. 
%
[R,C]=size(imn);

% downsampling: when the data is missed, we insert not a number (NaN) sign.
if scale_factor>1 && down_sample==1
    % INPAINTING
    valid_cols = [];
    nan_cols = [];
    emtyVector = zeros(R,1);
    emtyVector(emtyVector ==0) = 255;
    flag= 0;
    for jj = 1:C
         if (imn(:,jj) == emtyVector | flag < 0)
             nan_cols = [nan_cols jj];
             flag = flag + 1;
         else 
             valid_cols = [valid_cols jj];
             flag = flag -1;
         end
    end 
%     valid_cols=uint16(1:scale_factor:C);
%     nan_cols=1:C;
%     nan_cols(valid_cols)=[];
elseif scale_factor>1 && down_sample==0
    C=C*scale_factor;
    valid_cols=uint16(1:scale_factor:C);
    imn2=zeros(R,C);
    imn2(:,valid_cols)=imn;
    imn=imn2;% Now, size(imn)==size(imn2);
    clear imn2
    nan_cols=1:C;
    nan_cols(valid_cols)=[];
else
    % DENOISING 
    valid_cols = [];
    nan_cols = [];
    emtyVector = zeros(R,1);
    emtyVector(emtyVector ==0) = 255;
    flag= 0;
    for jj = 1:C
         if imn(:,jj) == emtyVector 
%         if imn(:,jj) == emtyVector | flag <= 0
             nan_cols = [nan_cols jj];
             flag = flag + 1;
         else
             valid_cols = [valid_cols jj];
             flag = flag -1;
         end
    end
%     nan_cols=[];
%     valid_cols=uint16(1:C);
end
imn(:,nan_cols)=nan;
strnumber=num2str(i);
outfname=['./datasets/ourDatasets/output/' strnumber];
imwrite(im2uint8(imn/255),[outfname '.tif'],'tif');
figure(1)
imagesc(imn)
colormap('gray')
%% extract patches from noisy image(imn) & its reduced noise version (imf)
% Xn & Xf
fprintf('extract patches from imn & imf ....\n');
% Xn=extract_patches_lex(imn,ps,step);
Xn=Get_patches_2_lex(imn,ps);

Xn=single(Xn);
nlabel=ones(1,size(Xn,2));
if ~isempty(im)
    imf=nan(size(imn));
    imf(:,valid_cols)=im(:,valid_cols);
else
    % ** Plug a Denoising Algorithm Here *************************************
%     psnlm=ps;
    imf=nan(size(imn));%---------
%    addpath('D:\rz');
    inp=imn(:,valid_cols);
%     Testp{1}=extract_patches_lex_col(inp,psnlm,step);
    Testp{1}=Get_patches_2_lex_col(inp,psnlm);
    [R2,C2]=size(inp);
    wei_arr=ones(1,size(Testp{1},2));
    Xf2= patchfiltering (Testp,wei_arr,psnlm(1), psnlm(2),...
        prod(psnlm),inp);
    imf(:,valid_cols)=insert_patches_Get_patches_2_lex_col(Xf2,R2,C2,psnlm);%-------
    clear Xf2 Testp inp wei_arr
    %~~~~~~~~~~~~~~~~~~
%     imf=nan(size(imn));
%     max_val=255;%max(max(imn(:,valid_cols)));
%     sigma=function_stdEst(im2uint8(imn(:,valid_cols)./max_val));
%     [NA, imf(:,valid_cols)] = BM3D(1, double(imn(:,valid_cols)/max_val), sigma);
%     imf(:,valid_cols)=imf(:,valid_cols)*max_val;
    %**********************************************************************
end
imf=single(imf);
% Xf=Get_patches_2_nofilter(imf,ps);
% Xf=extract_patches_lex(imf,ps,step);
Xf=Get_patches_2_lex(imf,ps);
Xf=single(Xf);
imwrite(im2uint8(imf/255),['denoisedimNLM600.tif'],'tif');
%% @@@@@ compute running time
time_start = cputime;
%%  remove mean of intensity 
% output: Xn
if NUM==1
    [Xn2,dc2]=remove_mean_inpainting_original(Xn);
    Xf2=remove_mean_inpainting_original(Xf);
else
    Xn2=remove_mean_inpainting_original(Xn);
    [Xf2,dc2]=remove_mean_inpainting_original(Xf);
end
%% find patterns of NaN
% Note that as we mentioned in the paper, we find similar patches for each
% patch that have similar pattern of missing values. In our experiments,
% this has better reconstruction quality.
r=ps(1);
c=ps(2);
if scale_factor>=1
%     Np=scale_factor;%number of patterns
%    p=ones(size(Xn,1),Np);
    count=0;
    jj=0;
    for ii=1:size(Xn,2)
        if sum(isnan(Xn(:,ii)))
            jj=jj+1;
            nan_patch(:,jj)=Xn(:,ii);
             count=count+1;
%             p_temp(:,count)=ones(64,1);
%             p_temp(:,count)=isnan(Xn(:,ii));
            p_nan(:,count)=ones(r*c,1);
            p_nan(:,count)=isnan(Xn(:,ii));
            p_temp(:,count)=ones(r*c,1);
            p_temp(:,count)=isnan(Xn(:,ii));
            for kk=1:count-1
                if p_nan(:,count)== isnan(p_temp(:,kk))
                     count=count-1;
                    break
                end
            end
            p_nan(p_nan(:,count)==1,count)=nan;
            p=ones(size(Xn,1),count);
            for ii=1:count
                p(:,ii)=p_nan(:,ii);
            end
            p_temp(p_temp(:,count)==1,count)=NaN;
        end
    end
end
p(:,count+1)=ones(r*c,1);
fprintf('# of nanpacthes: %d\n',jj);
fprintf('# of patterns: %d\n',size(p,2));
%% find similar patches for each patch
% find 'num' NL (nonlocal) patches for each patch & gather them into a big matrix (XX)
fprintf('find similar patches & put them into XX{k} ....\n');
[XX,nanpatches,patches,wei]=find_nl_for_inpainting(Xn2,Xf2,nlabel,1,NUM,NUM_neighbors,jj);
num_patches=size(Xn2,2);
clear Xn2 Xf2
%% OMP for inpainting
% load(dlfile)
fprintf('OMP for inpainting ....\n');
alpha=cell(1,1);
par_test{k}.num=NUM;
if scale_factor>=1
    par_test{k}.nan_patterns=p; 
end
% ** Sparse Coding method

alpha{k}=sparse_func(num_patches,XX{k},nanpatches{k},patches{k},D{k},par_test{k});
alpha{k}=sparse(alpha{k});

%% reconstruct patches in each cluster
fprintf('reconstruction ...\n');
Xnhat=zeros(size(Xn));
prob_out=single(zeros(2,size(Xn,2)));
%
id=nlabel==k;
N=sum(id);
% d=size(alpha{k},1);
h=80;
% Compute nonlocal weighted sparse representation (NWSR)
[Xnhat(:,id),prob_out(:,id)]=omp_mean_patches_wei_noDC(D{k},alpha{k},NUM,N,h,wei{k});
%% Add the removed mean of intensity to the patches.
% reconstructed patches are gathered into Xnhat. Thus we need to add mean
% of intensity (DC component) to these patches.
dc1=dc2.*prob_out(1,:)+dc2.*prob_out(2,:);
clear dc2;
Xnhat_dc=Xnhat+repmat(dc1,size(Xnhat,1),1);
% for ii=1:size(Xnhat_dc,2)
%     if any(isnan(Xnhat_dc(:,ii)))
%         error('nan patch idx is %d', ii);
%     end
% end
%% @@@@@
time_end = round(cputime-time_start);
%% Reconstruction of the whole image from the patches.
im_out=insert_patches_Get_patches_2_lex(Xnhat_dc,R,C,ps);
% toc
% clear Xn id
im_out=double(im_out);