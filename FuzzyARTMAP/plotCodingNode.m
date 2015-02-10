function out_nothing=plotCodingNode(input_point,dep_f_val,filesavedir,filename_prefix,sampleNumber,varargin)

check_cond=((sampleNumber>-320) && (sampleNumber<323));
%check_cond=1;

rect_created=[(min(input_point(1),1-input_point(3))) (min(input_point(2),1-input_point(4))) abs(1-input_point(3)-input_point(1))+.01 abs(1-input_point(4)-input_point(2))+.01];

if nargin==5
    out_nothing=rectangle('Position',rect_created,'EdgeColor',[0 0 0],'LineWidth',1.5);
    hold on
    tmp_handle=text((rect_created(1)+rect_created(3)/2),(rect_created(2)+rect_created(4)/2),'A_{dep}');
elseif nargin==6
    if varargin{1}==2
        out_nothing=rectangle('Position',rect_created,'EdgeColor',[.6 .6 1],'LineWidth',1.5);
        hold on

    else
        out_nothing=rectangle('Position',rect_created,'EdgeColor',[1 .6 .6],'LineWidth',1.5);
        hold on

    end
elseif nargin==7
    out_nothing=rectangle('Position',rect_created,'EdgeColor',[0 0 0],'LineWidth',1.5);
    hold on
    tmp_handle=text((rect_created(1)+rect_created(3)/2),(rect_created(2)+rect_created(4)/2),'w_{current}') ;

end

if check_cond
    if nargin==7
        xlabel([regexprep(get(get(gca,'xlabel'),'string'),'curr|updtd','') 'curr R_j: [' regexprep(num2str(round(input_point*100)'/100),'\s+',' '),'] |R_j| ' num2str(round((2-sum(input_point))*1000)/1000)])
    elseif nargin==6
        xlabel([regexprep(get(get(gca,'xlabel'),'string'),'curr|updtd','') ' updtd R_j: [' regexprep(num2str(round(input_point*100)'/100),'\s+',' '),'] |R_j| ' num2str(round((2-sum(input_point))*1000)/1000)])
    elseif nargin==5
        xlabel([regexprep(get(get(gca,'xlabel'),'string'),'curr|updtd','') ...
            ' curr f_{dep}: [' regexprep(num2str(round(dep_f_val*100)'/100),'\s+',' '),']'...
            ' curr A_{dep}: [' regexprep(num2str(round(input_point*100)'/100),'\s+',' '),']'])
    end


    axis equal
    axis tight
    axis([-.1 1.1 -.1 1.1]);

    ginput(1);


end

if exist('tmp_handle','var')
    delete(tmp_handle)
end
if nargin==7
    delete(out_nothing)

    if varargin{1}==2
        out_nothing=rectangle('Position',rect_created,'EdgeColor',[.6 .6 1],'LineWidth',1.5);
        hold on
    else
        out_nothing=rectangle('Position',rect_created,'EdgeColor',[1 .6 .6],'LineWidth',1.5);
        hold on
    end

end

% else
%     out_nothing='';
% end
