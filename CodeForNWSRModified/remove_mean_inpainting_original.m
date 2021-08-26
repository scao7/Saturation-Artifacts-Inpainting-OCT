% this function removes mean of intensity values. If inputs have missing
% values, this function neglects them and comptutes mean of intensity based
% on observed values.
%
% select_dimension=1, the signals are stored in columns. 2, the signals are
% stored in rows
% 
function [Xn_noDC,mean_columns]=remove_mean_inpainting_original(Xn,select_dimension)
if ~exist('select_dimension','var') || isempty(select_dimension)
    select_dimension=1;
end
if select_dimension==1
    num_signals=size(Xn,2);
    for ii=1:num_signals
        dim=size(Xn,1);
        num_nans=sum(isnan(Xn(:,ii)));
        num_valid=dim-num_nans;
        if num_nans~=0
            Xn2=Xn(:,ii);
            Xn2=reshape(Xn2(~isnan(Xn2)),num_valid,1);
        else
            Xn2=Xn(:,ii);
        end
        mean_columns(1,ii)=mean(Xn2);
        t=~isnan(Xn(:,ii));
        Xn_temp=nan(size(Xn(:,ii)));
        Xn_temp(t)=Xn2-repmat(mean_columns(1,ii),num_valid,1);
        Xn_noDC(:,ii)=Xn_temp;
    end
%     dim=size(Xn,1);
%     num_nans=sum(isnan(Xn(:,1)));
%     num_valid=dim-num_nans;% number of valid elements 
%     mean_columns=mean(Xn2);
%     t=~isnan(Xn);
%     Xn_noDC=nan(size(Xn));
%     Xn_noDC(t)=Xn2-repmat(mean_columns,num_valid,1);% size(Xn2,1),1
else
    % signals are stored in rows of the matrix Xn
    num_signals=size(Xn,1);
    dim=size(Xn,2);
    num_nans=sum(isnan(Xn(1,:)));
    num_valid=dim-num_nans;% number of valid elements
    if num_nans~=0
        Xn2=Xn;
        Xn2=reshape(Xn2(~isnan(Xn2)),num_signals,num_valid);
    else
        Xn2=Xn;
    end
    mean_columns=mean(Xn2,2);
    t=~isnan(Xn);
    Xn_noDC=nan(size(Xn));
    Xn_noDC(t)=Xn2-repmat(mean_columns,1,num_valid);
end