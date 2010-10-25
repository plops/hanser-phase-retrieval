% opposite of decimate:
% repeat each pixel of the input image zoom times
% (like resample but only works with integral values zoom
% and doesn't do interpolation)

function big=decimate(in,zoom)
big_m=zeros(zoom*size(in));
in_m=double(in);
[row col]=size(big_m);
% let matlabs jit compiler work its magic
for j=1:col
    for i=1:row
        big_m(i,j)=in_m(ceil(i/zoom),ceil(j/zoom));
    end
end
big=dip_image(big_m);
end