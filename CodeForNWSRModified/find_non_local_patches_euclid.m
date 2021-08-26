% find non local patches
%   it just finds nearest patches with the same label that are close to
%   X(:,idx)
%   Assume: row-wise patch extraction method e.g. EXTRACT_PATCHES_LEX.M
% EXAMPLE
%     addpath('E:/THESIS/Implements/Pedagogical/ML/');
% Ashkan
function [nlidx,euc_out]=find_non_local_patches_euclid(X,idx,label,num_neighbors,num)
nlidx=find_non_local_patches_fast(idx,label,num_neighbors);
% nlidx=find_non_local_patches_fast_clustering(idx,label,num_neighbors,2*num_neighbors);% L>=2*num_neighbors
%
C=X(:,idx)';
XX=X(:,nlidx)';
%
% r=ps(1);
% c=ps(2);
f=[];
if any(isnan(C))
    % select the patches with similar NaN patterns to X(:,idx)
    pat=isnan(C);% NaN patterns in the patch
    mask=repmat(pat,size(XX,1),1);
    XX(mask==1)=nan;
    loc=sum(isnan(XX)==repmat(pat,size(XX,1),1),2);
%     id=loc~=0;
    id=loc==64;
    XX2=XX(id,:);
    % remove NaNs
    C=remove_nans(C,2);
    XX2=remove_nans(XX2,2);
    f=find(id);
    %---------------
% % % %     C=remove_nans(C,2);
% % % %     XX2=remove_nans(XX,2);
% % % %     f=1:size(XX,1);
else 
    XX2=XX;
end
%
param.sort=1;
[cls_idx,euc]=euclid_classifier(XX2,C,param,idx);
% [cls_idx,euc]=mahal_classifier(XX2,C,param);

% [cls_idx,euc]=mahal_classifier_2_mex(XX2,C);% fastest, 35 seconds!

% [cls_idx,euc]=mahal_classifier_2_mex(XX2,mean(XX2,1));
% [cls_idx,euc]=mahal_classifier_2(XX2,C); % faster
%
if isempty(f)
    nlidx=nlidx(cls_idx(1:num));
%     nlidx=[1;2;3;4;5;6];
    euc_out=euc(cls_idx(1:num));
%     euc_out=repmat(euc_out,6,1);
else
    nlidx=nlidx( f(cls_idx(1:num)));% X(:,nlidx(f(cls_idx)))'
    if nargout>1
        euc_out=euc(cls_idx(1:num));
    end
end

for ii=1:num
    if isnan(euc_out(ii))
        error('nan euc. The idx is %d\n', idx);
    end
end