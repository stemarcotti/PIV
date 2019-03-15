function plot_data = generate_plot_normalized(x_field, y_field, scale)

[y_size, x_size] = size(x_field);
i = 0;
for x = 1:scale:x_size
  for y = y_size:-1*scale:1
    if (isnan(x_field(y, x)) == 0)
      r = sqrt(x_field(y, x)^2 + y_field(y, x)^2);
      i = i + 1;
      plot_data(i, :) = [x, y, x_field(y, x) / r, y_field(y, x) / r];
    end
  end
end