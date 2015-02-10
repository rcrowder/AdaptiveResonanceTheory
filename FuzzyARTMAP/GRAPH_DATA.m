function A = GRAPH_DATA(train_data,train_out_desired,train_out_obtained,varargin)

if nargin==3
DrawPoint='.';
else
    DrawPoint=varargin{1};
end


color_red_right=[1 .7 .7];
color_red_wrong=[1 0 0];
color_blue_right=[.7 .7 1];
color_blue_wrong=[0 0 1];

color_red_right=[1 .4 .4];
color_red_wrong=[1 0 0];
color_blue_right=[.4 .4 1];
color_blue_wrong=[0 0 1];
% 
% color_red_right=[1 0 0];
% color_red_wrong=[1 0 0];
% color_blue_right=[0 0 1];
% color_blue_wrong=[0 0 1];


for i=1:length(train_data)
    if train_out_obtained(i)==1
        if train_out_desired(i)==1
            plot(train_data(1,i),train_data(2,i),DrawPoint,'Color',color_red_right,'markersize',14);
            hold on;
        else
            plot(train_data(1,i),train_data(2,i),DrawPoint,'Color',color_red_wrong,'markersize',14);
            hold on;
        end
    else
        if train_out_desired(i)==2
            plot(train_data(1,i),train_data(2,i),DrawPoint,'Color',color_blue_right,'markersize',14);
            hold on;
        else
            plot(train_data(1,i),train_data(2,i),DrawPoint,'Color',color_blue_wrong,'markersize',14);
            hold on;
        end
    end
    %Use this when showing ordinals in isolation
    %text(train_data(1,i)-.055,train_data(2,i),num2str(i),'fontname','arial','fontsize',13);
    
    %SUBPLOT Use this when stepping through ART learning
    %text(train_data(1,i)-.035,train_data(2,i)+.045,num2str(i),'fontname','
    %arial','fontsize',13,'color',[.5 .5 .5]);
    
    %PLOT Use this when stepping through ART learning
    %text(train_data(1,i)-.055,train_data(2,i)+.045,num2str(i),'fontname','arial','fontsize',13,'color',[.5 .5 .5]);
    
end

