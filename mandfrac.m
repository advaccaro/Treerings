% Mandelbrot set fractal generator
% Ported by Carl M (http://blog.eop.org.uk)
% Original by Kin-Hing Lo (http://cow.neondragon.net)

clear all;
clc;
clf;

img_height = 1000; % set image size
img_width = img_height;

x_min = -2.1;
x_max = 0.6;
y_min = -1.2;
y_max = 1.2;

fraimg = zeros(img_width,img_height);

i_max = 18; % number of iterations to do/colours - higher is better
starttime = clock;

for j=1:img_height

    for i=1:img_width

        x = (x_min+(i*((x_max - x_min) / (img_width - 1))));
        y = (y_min+(j*((y_max - y_min) / (img_height - 1))));
        
        iteration = 1;
        z = [0,0];
        x2 = 0;
        y2 = 0;

        while ( (iteration <= i_max) && ((x2)+(y2) <= 4))
           
            k = z;

            z(1) = (x2) - (y2);
            z(2) = 2*k(1)*k(2);

            z(1) = z(1) + x;
            z(2) = z(2) + y;

            x2 = (z(1))*(z(1));
            y2 = (z(2))*(z(2));

            iteration = iteration + 1;           
        end
        if iteration == i_max
            fraimg(i,j) = i_max + 1; % black inside fractal
        else
            fraimg(i,j) = iteration;
        end        
    end
end
colours = colormap(hot(i_max - 1));
colours(i_max,1:3) = 0; % create black
colormap(colours);
endtime = clock;
image(fraimg);
axis off equal image;
text(10,10,num2str(etime(endtime,starttime)),'BackgroundColor',[.7 .9 .7]);
