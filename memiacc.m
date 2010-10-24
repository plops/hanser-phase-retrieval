%%
vars;
acc=newim(s)*1i;
for phi=-alpha:alpha/10:alpha
    for theta=0:alpha/10:alpha
        acc=acc+memi(phi,theta,rr(s)<20,rr(s)<30);
    end
end
