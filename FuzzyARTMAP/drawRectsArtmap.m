for i=1:size(c.w,2)
    rect_created=[(min(c.w(1,i),1-c.w(3,i))) (min(c.w(2,i),1-c.w(4,i))) abs(1-c.w(3,i)-c.w(1,i))+.01 abs(1-c.w(4,i)-c.w(2,i))+.01];
     
    %IF RESCALING
    %rect_created=[(min((9*c.w(1,i)-4),1-(9*c.w(3,i)-4))) (min((9*c.w(2,i)-4),1-(9*c.w(4,i)-4))) abs(1-(9*c.w(3,i)-4)-(9*c.w(1,i)-4))+.01 abs(1-(9*c.w(4,i)-4)-(9*c.w(2,i)-4))+.01];

    if c.W(i,1)==0
        rectangle('Position',rect_created,'EdgeColor',[0 0 1],'LineWidth',1.5)
        hold on
        %text((rect_created(1)+rect_created(3)/2),(rect_created(2)+rect_created(4)/2),'\beta')
    else
        rectangle('Position',rect_created,'EdgeColor',[1 0 0],'LineWidth',1.5)
        hold on
        %text((rect_created(1)+rect_created(3)/2),(rect_created(2)+rect_created(4)/2),'\alpha')
    end
end

