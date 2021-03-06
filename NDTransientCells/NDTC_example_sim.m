function [X,Z,B]=NDTC_example_sim(Lstr,Ldur,B1,K2)

X=[0];
Z=[1];

dt=0.001; % 1 ms
Tb=0.1;

x_t1=0;
z_t1=1;
for t=2:500
    if t<=Ldur
        I_t1=Lstr;
    else
        I_t1=0;
    end
    [x_t2,z_t2,B_t2]=NDTransientCells(x_t1,z_t1,I_t1,B1,K2,Tb);
    X=[X x_t2];
    Z=[Z z_t2];
    x_t1=x_t2;
    z_t1=z_t2;
end

B=max(X.*Z-Tb,0);

% figure
% % axis([1 1000 0 1])
% plot(X)
% xlabel('Time (ms)','Fontsize',12)
% ylabel('Activity','Fontsize',12)
% title(['Change-sensitive receptor'],'Fontsize',12)
% box off
% 
% figure
% % axis([1 1000 0 1])
% plot(Z,'g')
% xlabel('Time (ms)','Fontsize',12)
% ylabel('Activity','Fontsize',12)
% title(['Habituative transmitter'],'Fontsize',12)
% box off
% 
% figure
% % axis([1 1000 0 1])
% plot(B,'r')
% xlabel('Time (ms)','Fontsize',12)
% ylabel('Activity','Fontsize',12)
% title(['Non-directional transient cell'],'Fontsize',12)
% box off

return