function colormap_2(low, up)
color_low  = low;
color_up = up;

N = 100;
for ii =1:3
    colorAll(:, ii) = linspace(color_low(ii), color_up(ii), N);
end
colormap(colorAll);

end