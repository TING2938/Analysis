function colormap_3(low, mid, up)
color_low  = low;
color_mid = mid;
color_up = up;

N = 100;
for ii =1:3
    colorAll_low(:, ii) = linspace(color_low(ii), color_mid(ii), N);
    colorAll_up(:, ii) = linspace(color_mid(ii), color_up(ii), N);
end
colorAll = [colorAll_low; colorAll_up];
colormap(colorAll);

end