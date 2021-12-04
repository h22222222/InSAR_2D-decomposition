% 二维形变提取
% 原理：最小二乘法
% 前提：一定要有升轨和降轨的数据，若只有单轨会出现秩亏现象，造成解不稳定
% revised by yuqinghe 2021/12/1
%% 得到三方向在LOS向的投影系数
% ALOS2 升轨
% Mean heading angle: -10.4586 degrees
% Mean incidence angle: 36.3027 degrees
set(0,'defaultfigurecolor','w');
ha_AT=-10.4586/180*pi;
ia_AT=36.3027/180*pi;
% 降轨
% Mean heading angle:-164.583degrees
% Mean incidence angle: 33.9736 degrees
ha_DT=-164.583/180*pi;
ia_DT= 33.9736/180*pi;
tem=3*pi/2;
% 升轨
prN_AT=-sin(ia_AT)*cos(ha_AT-tem);
prE_AT=-sin(ia_AT)*sin(ha_AT-tem);
prU_AT=cos(ia_AT);
% 降轨
prN_DT=-sin(ia_DT)*cos(ha_DT-tem);
prE_DT=-sin(ia_DT)*sin(ha_DT-tem);
prU_DT=cos(ia_DT);
% LOS vector [ENU]
los_AT = [prN_AT prE_AT prU_AT];
los_DT = [prN_DT prE_DT prU_DT];
% los = [0.3 -0.15 .9];%ENVISAT
% los = [0.5487 -0.1000 0.8300]; %terraSAR descending
% los_az = [-0.1793 -0.9838 0];%terraSAR descending in azimuth
%los = [ -0.4150 -0.0791 0.9064]; ascending
% Poissions ratio
nu = 0.25;

% inputFile='mengu2021DT04_alldef.mat';
atinputFile=load('mengu2021ALOS_downsampled.txt');
% 计算需要，去掉形变为0/特定值 的行
% [judge,~]=find(atinputFile==0);
% if ~isempty(judge)
% 
%        atinputFile(unique(judge),:)=[];
% 
% end
    

at_lon=atinputFile(:,2);
at_lat=atinputFile(:,1);
at_obs=atinputFile(:,3);
at_length=length(at_lon);
dtinputFile=load('menggu_dt04_downsampled.txt');
% [judge,~]=find(dtinputFile==0);
% if ~isempty(judge)
% 
%        dtinputFile(unique(judge),:)=[];
% 
% end
dt_lon=dtinputFile(:,2);
dt_lat=dtinputFile(:,1);
dt_obs=dtinputFile(:,3);
dt_length=length(dt_lon);

%% 将升降轨点位划归到相同点上，如果不存在大面积失相干的区域可采用插值得到相同点，反之可认为在某个经纬度差异下，相近的点位是一致的。
%当点位的差异不超过0.001度时，认为是同一个点，可根据最后的数据点个数做调整  at的范围小一点
result_at=[];result_dt=[];
for i=1:at_length
    for j=1:dt_length
        if (abs(at_lon(i)-dt_lon(j))<0.001)&(abs(at_lat(i)-dt_lat(j))<0.001)
            result1=[at_lon(i) at_lat(i) at_obs(i)];
            result2=[at_lon(i) at_lat(i) dt_obs(j)];
            result_at=[result_at;result1];
            result_dt=[result_dt;result2];
        end
    end
end

%绘图
figure;
scatter(result_at(:,1),result_at(:,2),95,result_at(:,3),'.');
colormap jet;
colorbar;
grid on;
title('AT145 def');
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')

figure;
scatter(result_dt(:,1),result_dt(:,2),95,result_dt(:,3),'.');
colormap jet;
colorbar;
grid on;
title('DT04 def');
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')

%% 最小二乘分解
dlosat=result_at(:,3);
dlosdt=result_dt(:,3);
n=length(dlosdt);
% D=[dlos1;dlos2];

X_PARA=[];
A_1=[los_AT(2) los_AT(3);los_DT(2) los_DT(3)];
for i=1:n
D_1=[dlosat(i,:);dlosdt(i,:)];
% X_PARA=inv(A)*B_NEW;
% x = lsqlin(C,d,A,b) 在满足 A*x ≤ b 的情况下基于最小二乘思想求解线性方程组 C*x = d。

X_PARA_1=lsqlin(A_1,D_1,[],[]);
X_PARA=[X_PARA;X_PARA_1];
end

d_e=[];
d_v=[];
for i=1:n
d_e1=X_PARA(2*i-1,1); 
d_v1=X_PARA(2*i,1); 
d_e=[d_e;d_e1];d_v=[d_v;d_v1];
end

%将结果绘图
figure;
scatter(result_dt(:,1),result_dt(:,2),125,d_e,'.');
colormap jet;
colorbar;
grid on;
title('def in EW direction','FontSize',15);
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
set(get(colorbar,'Title'),'string','def(cm)');

figure;
scatter(result_dt(:,1),result_dt(:,2),125,d_v,'.');
colormap jet;
colorbar;
grid on;
title('def in vertical direction','FontSize',15);
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
set(get(colorbar,'Title'),'string','def(cm)');

%% 输出到文件
grdwrite2(result_dt(:,1),result_dt(:,2),d_v,'vertical_def.grd');
grdwrite2(result_dt(:,1),result_dt(:,2),d_e,'EW_def.grd');

fid=fopen('vertical_def.out','w');
[m,n]=size(result_dt);
for i=1:m
       
        fprintf(fid,'%f\t%f\t%f\n',result_dt(i,1),result_dt(i,2),d_v(i,:));

end
fclose(fid);
        
fid=fopen('EW_def.out','w');
[m,n]=size(result_dt);
for i=1:m
       
        fprintf(fid,'%f\t%f\t%f\n',result_dt(i,1),result_dt(i,2),d_e(i,:));

end
fclose(fid);

%% 绘制跨过断层的剖线进行分析
%东西向形变场
figure;
scatter(result_dt(:,1),result_dt(:,2),125,d_e,'.');
colormap jet;
colorbar;
grid on;
title('def in EW direction','FontSize',15);
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
set(get(colorbar,'Title'),'string','def(cm)');
hold on;
% plot([x1(i),x2(i)],[y1(i),y2(i)],'r-')

% plot([99.5,101],[51.25,51.25],'.-k','LineWidth',0.05);
% plot([99.5,101],[51.17,51.17],'.-k','LineWidth',0.05);

plot([100.1,100.5],[51.55,51.55],'.-k','LineWidth',0.05);
plot([100.1,100.5],[51.48,51.48],'.-k','LineWidth',0.05);
plot([100.1,100.5],[51.4,51.4],'.-k','LineWidth',0.05);
plot([100.1,100.5],[51.34,51.34],'.-k','LineWidth',0.05);
hold on;
text(100.5,51.55,'A''');text(100.1,51.55,'A');
text(100.5,51.48,'B''');text(100.1,51.48,'B');
text(100.5,51.4,'C''');text(100.1,51.4,'C');
text(100.5,51.34,'D''');text(100.1,51.34,'D');
% h = colorbar;
% h.Label.String = 'deformation (m)';
% set(h,'ytick',[-0.06:0.02:0.065]);
% axis image;axis xy;
% set(gca,'XLim',[76.8 77.6]);%X轴的数据显示范围
hold on;

linedeA=improfile();%指针对imagesc的figure有效，不行的时候可以采用经纬度取值
save('./defprofile/linedeA.mat','linedeA');

result_de_A=[];result_dv_A=[];
result_de_B=[];result_dv_B=[];
result_de_C=[];result_dv_C=[];
result_de_D=[];result_dv_D=[];
delength=length(d_e);dvlength=length(d_v);
% [100.1,100.5],[51.55,51.55],
% [100.1,100.5],[51.48,51.48]
% [100.1,100.5],[51.4,51.4]
% [100.1,100.5],[51.34,51.34]
for i=1:delength  %abs(result_dt(i,1)-100.1)<0.008)
        if abs(result_dt(i,2)-51.34)<0.001
            result1=[result_dt(i,1) result_dt(i,2) d_e(i)];
%             result2=[at_lon(i) at_lat(i) dt_obs(j)];
             result_de_D=[result_de_D;result1];
%             result_dt=[result_dt;result2];
        end
end



%垂直向形变场
figure;
scatter(result_dt(:,1),result_dt(:,2),125,d_v,'.');
colormap jet;
colorbar;
grid on;
title('def in vertical direction','FontSize',15);
xlabel('Longitude (degrees)')
ylabel('Latitude (degrees)')
set(get(colorbar,'Title'),'string','def(cm)');
hold on;
% plot([x1(i),x2(i)],[y1(i),y2(i)],'r-')

% plot([99.5,101],[51.25,51.25],'.-k','LineWidth',0.05);
% plot([99.5,101],[51.17,51.17],'.-k','LineWidth',0.05);

plot([100.1,100.5],[51.55,51.55],'.-k','LineWidth',0.05);
plot([100.1,100.5],[51.48,51.48],'.-k','LineWidth',0.05);
plot([100.1,100.5],[51.4,51.4],'.-k','LineWidth',0.05);
plot([100.1,100.5],[51.34,51.34],'.-k','LineWidth',0.05);

text(100.5,51.55,'A''');text(100.1,51.55,'A');
text(100.5,51.48,'B''');text(100.1,51.48,'B');
text(100.5,51.4,'C''');text(100.1,51.4,'C');
text(100.5,51.34,'D''');text(100.1,51.34,'D');
% text(76.82,40.25,'(a)','fontsize',18);
% set(gca,'YDir','normal'); 
colormap jet;
colorbar;
h = colorbar;
h.Label.String = 'deformation (m)';
% set(h,'ytick',[-0.06:0.02:0.065]);
% axis image;axis xy;
% set(gca,'XLim',[76.8 77.6]);%X轴的数据显示范围

%绘制剖线结果
figure;
subplot(2,2,1);
% yyaxis left;
% plot(lineAT772_new*100','.','LineWidth',0.05);
% hold on;
plot(result_de_A(:,3),'-','LineWidth',0.5);
hold on;
plot(result_de_A(:,3),'.','LineWidth',0.05);
ylabel('Def. (cm)');
title('Deformation profile along AA*');
xlabel('Along the profile');
% ylabel('Elev. (m)');
grid on;
subplot(2,2,2);
% yyaxis left;
% plot(lineAT772_new*100','.','LineWidth',0.05);
% hold on;
plot(result_de_B(:,3),'-','LineWidth',0.5);
hold on;
plot(result_de_B(:,3),'.','LineWidth',0.05);
ylabel('Def. (cm)');
title('Deformation profile along BB*');
xlabel('Along the profile');
% ylabel('Elev. (m)');
grid on;
subplot(2,2,3);
% yyaxis left;
% plot(lineAT772_new*100','.','LineWidth',0.05);
% hold on;
plot(result_de_C(:,3),'-','LineWidth',0.5);
hold on;
plot(result_de_C(:,3),'.','LineWidth',0.05);
ylabel('Def. (cm)');
title('Deformation profile along CC*');
xlabel('Along the profile');
% ylabel('Elev. (m)');
grid on;
subplot(2,2,4);
% yyaxis left;
% plot(lineAT772_new*100','.','LineWidth',0.05);
% hold on;
plot(result_de_D(:,3),'-','LineWidth',0.5);
hold on;
plot(result_de_D(:,3),'.','LineWidth',0.05);
ylabel('Def. (cm)');
title('Deformation profile along DD*');
xlabel('Along the profile');
% ylabel('Elev. (m)');
grid on;