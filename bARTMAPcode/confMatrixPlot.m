%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of an ARTMAP variant,
% as described in:
% Gaddam, C. S. (2007).
% Feature Selection via transmitter depletion in ARTMAP. Online Document, xx(x) xxx-xxx.
% Boston, MA: Boston University.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Sai Chaitanya Gaddam (August 2007-08)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Usage: [conf_actual,conf_normalized]=confMatrix(actual,predicted|,[class_labels]|,dataset_title)
%Strings in class_labels will be used for axes labels if provided
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A,A_norm,avg_acc]=confMatrixPlot(actual,predicted,c,varargin)

A=zeros(size(max(length(actual),length(unique(predicted)))));
A_norm=zeros(size(max(length(actual),length(unique(predicted)))));

if (nargin==5) || ((nargin==4) && ischar(varargin{1}))
    dataset_title=[varargin{nargin-3}, ': '];
else
    dataset_title='';
end

fontsize_val=12;
fontsize_val_offdiag=12;
fontsize_val_labels=12;
label_rotation=45;


%For Genome Data
% fontsize_val=4;
% fontsize_val_offdiag=4;
% fontsize_val_labels=4;


%actual_unique=unique(actual);
actual_unique=unique([c.W*[1:1:size(c.W,2)]' ; unique(actual)]);

for i=1:length(actual_unique)
    for j=1:length(actual_unique)
        A(i,j)=sum((actual(:)==actual_unique(i)).*(predicted(:)==actual_unique(j)));
        A_norm(i,j)=sum((actual(:)==actual_unique(i)).*(predicted(:)==actual_unique(j)))/sum((actual(:)==actual_unique(i)));
    end
end



figure(1)
imagesc(A*0)
%colormap([[1:-.01:.7]' [1:-.01:.7]' [1:-.01:.7]'] );
colormap([[1:-.01:1]' [1:-.01:1]' [1:-.01:1]'] );

class_labels=cell(length(actual_unique),1);
if ((nargin==3) || ((nargin==4) && ischar(varargin{1})))
    for i=1:length(actual_unique)
        class_labels{i}=num2str(actual_unique(i)-1);
    end
else
    class_labels=varargin{1}(actual_unique);
end

axis equal;
axis tight;

%set(gcf,'position',[1 -419 1600 946])
set(gca,'xticklabel','');
set(gca,'yticklabel','');
for i=1:length(actual_unique)
    for j=1:length(actual_unique)
        if (round(A(j,i))>0)
            if  (i~=j)
                h_curr=text(i,j,num2str(A(j,i)),'fontsize',fontsize_val_offdiag,'color',0*[0.4 0.4 0.4],'horizontalalignment','right');
            else
                h_curr=text(i,j,num2str(A(j,i)),'fontsize',fontsize_val,'fontweight','bold','color',[0 0 0],'horizontalalignment','right');
                h_extent=get(h_curr,'extent');
                text(h_extent(1)+h_extent(3),j,[' /' num2str(sum(A(j,:)))],'fontsize',fontsize_val,'color',[0 0 0],'horizontalalignment','left');
            end
        end
        if (j==1)
            if ~isnan(A_norm(i,i))
                percent_corr=num2str(A_norm(i,i)*100,3.2);
            else
                percent_corr='no test pts.';
            end

            text(length(actual_unique)+1,i,percent_corr,'fontsize',fontsize_val,'color',[0.8 0 0],'horizontalalignment','center');
            text(0.4,i,class_labels{i},'fontsize',fontsize_val_labels,'color',[0 0 0],'horizontalalignment','right','rotation',label_rotation);
        end
        if (i==1)
            %text(length(actual_unique)+1,i,num2str(sum(A(i,:))),'fontsize',fontsize_val,'color',[.8 0 0],'horizontalalignment','center');
            text(j,length(actual_unique)+.6,class_labels{j},'fontsize',fontsize_val_labels,'color',[0 0 0],'horizontalalignment','right','rotation',label_rotation);
        end

    end
end



text(length(actual_unique)+.75,0.5,'% Correct:' ,'fontsize',fontsize_val,'color',[.2 .2 .2],'horizontalalignment','left');
text(length(actual_unique)+.75,length(actual_unique)+.5,'Mean Accuracy:' ,'fontsize',fontsize_val,'color',[.2 .2 .2],'horizontalalignment','left');

A_norm_not_nan=diag(A_norm);
A_norm_not_nan=A_norm_not_nan(~isnan(A_norm_not_nan));
avg_acc=sum(A_norm_not_nan)*100/length(A_norm_not_nan);
avg_acc_str=num2str(avg_acc,3.2);
text(length(actual_unique)+ 1,length(actual_unique)+.75,[ avg_acc_str ' %'] ,'fontsize',fontsize_val,'color',[0.8 0 0],'horizontalalignment','center');

%h_title=title(dataset_title,'fontsize',fontsize_val);
%set(h_title,'position',get(h_title,'position'));
%h_ylab=ylabel('Actual Class','fontsize',fontsize_val_labels,'color',[.4 .4 1]);
%set(h_ylab,'position',get(h_ylab,'position')+[-0.5 0 0]);
%h_xlab=xlabel('Predicted Class','fontsize',fontsize_val_labels,'color',[.4 .4 1]);
%set(h_xlab,'position',get(h_xlab,'position')+[0 0.25 0]);
%x_pos=get(h_xlab,'position')+[0 0.25 0];
%text(x_pos(1),x_pos(2)+.2,['No. committed nodes: ' num2str(c.C) ' No. training pts: ' num2str(length(c.NodeList)) ],'fontsize',fontsize_val,'color',[0 0 0],'horizontalalignment','center');
% field_names=fieldnames(c);
% for i=1:length(field_names)
%
%     if length(c.(field_names{i}))<2
%     text(x_pos(1),x_pos(2)+.20*i,[field_names{i} ': ' num2str(c.(field_names{i}))],'fontsize',fontsize_val,'color',[0 0 0],'horizontalalignment','center');
%     else
%     text(x_pos(1),x_pos(2)+.20*i,[field_names{i} 'size : [' num2str(size(c.(field_names{i}))) ']'],'fontsize',fontsize_val,'color',[0 0 0],'horizontalalignment','center');
%     end
% end
set(gca,'ytick',[1:1:size(c.W,2)]+.5)
set(gca,'xtick',[1:1:size(c.W,2)]+.5)
set(gca,'linewidth',1.5)
grid

