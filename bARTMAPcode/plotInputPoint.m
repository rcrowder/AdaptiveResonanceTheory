function output_handle = plotInputPoint(input_point,filesavedir,filename_prefix,sampleNumber,input_class)
DrawPoint='.';
check_cond=((sampleNumber>-320) && (sampleNumber<323));
%check_cond=1;
if check_cond
    if input_class==1
        output_handle= plot(input_point(1),input_point(2),DrawPoint,'Color',[1 0 0]);
        tmp_handle=text(input_point(1),input_point(2),['input ',num2str(sampleNumber) ],'fontsize',8);
        hold on;
    else
        output_handle= plot(input_point(1),input_point(2),DrawPoint,'Color',[0 0 1]);
        tmp_handle=text(input_point(1),input_point(2),['input ',num2str(sampleNumber) ],'fontsize',8);
        hold on;
    end
    axis equal
    axis tight
    axis([-.1 1.1 -.1 1.1]);

    %     if ((sampleNumber>470) && (sampleNumber<478))
    %         input('Press Enter');
    %     end


    if check_cond
        if ~exist('tmp_a','var')
            [tmp_a,tmp_b,tmp_c]=mkdir('C:\Program Files\MATLAB\R2007b\work\',filesavedir);
        end
        filesavedir=['C:\Program Files\MATLAB\R2007b\work\' filesavedir];
        load k_stepTrack k_step
        saveas(gcf,[filesavedir filename_prefix '_' num2str(k_step)],'emf')
        h_temp=get(gca,'Children');
        hgsave(h_temp,[filesavedir filename_prefix '_fig_' num2str(k_step)])
        k_step=k_step+1;
        save k_stepTrack k_step
    end

    delete(tmp_handle)
end
figure(1)

% else
%
%     output_handle ='';
% end


