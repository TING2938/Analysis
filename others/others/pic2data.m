%% 提取图片坐标；
%%
path=uigetfile('.png');
rawData=imread(path);
rawData=rawData(:,:,1);
[b1,b2]=find(rawData<255); % 设置阈值
b1=max(b1)+10-b1;
scatter(b2,b1,'.');
[X,Y]=ginput(4);
%% 输入图中坐标；
x1=        0    ;
x2=          8.98403   ;
y1=            0;
y2=            4.49203  ;
%%
m=b1;
n=b2;
for i=1:length(b1)
    m(i)=x1+(b2(i)-X(1))*(x2-x1)/(X(2)-X(1));
    n(i)=y1+(b1(i)-Y(3))*(y2-y1)/(Y(4)-Y(3));
end
fig=scatter(m,n,'.','k');