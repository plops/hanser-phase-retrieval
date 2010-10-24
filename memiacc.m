%%
tic
vars;
acc=newim([s 600])*1.0;
ma=1/(2*D);
i=0;
for kx=-ma:ma/10:ma
    for ky=-ma:ma/10:ma
        i=i+1;
        acc(:,:,i)=squeeze(acc(:,:,i-1))+...
            abs(memi(kx,ky,rr(s)<20,rr(s)<30)).^2;
    end
end
toc
