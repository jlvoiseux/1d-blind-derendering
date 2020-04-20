using Images
using TestImages

function generate_conv_image()
    f = testimage("mandrill");
    y = imfilter(f, reflect(Kernel.gaussian(3)));
    imshow(y)
end

generate_conv_image()
