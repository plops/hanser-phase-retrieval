%% incoherent addition of images in focal plane
tic
s=[128 128]
acc=newim([s 600])*1.0;
ma=1/2;
i=0;
for kx=-ma:ma/10:ma
    for ky=-ma:ma/10:ma
        i=i+1;
        acc(:,:,i)=squeeze(acc(:,:,i-1))+...
            abs(memi(kx,ky,rr(s)<10,rr(s)<30)).^2;
    end
end
toc


