function res = pic2data(fileName, x1, x2, y1, y2)
%% get data from picture; 
    rawData = imread(fileName);
    thres = graythresh(rawData);
    rawData = im2bw(rawData, thres);
    [b1, b2] = find(rawData == 0);
    b1 = max(b1) + 10 - b1;
    scatter(b2, b1, '.k');
    [X,Y] = ginput(4);
%% 
    m = b1;
    n = b2;
    for i=1:length(b1)
        m(i)=x1 + (b2(i) - X(1)) * (x2 - x1) / (X(2) - X(1));
        n(i)=y1 + (b1(i) - Y(3)) * (y2 - y1) / (Y(4) - Y(3));
    end
    res.x = m;
    res.y = n;
    scatter(m, n, '.k');

end