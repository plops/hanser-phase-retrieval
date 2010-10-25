% calculate in-focus image for one coherent beam
% the direction of the beam specified direction is
% given by kx and ky both can assume values in [-.5 .. .5]
% example: vars; memi(0,-.5,rr(s)<20,rr(s)<20)

function u5 = memi(kx,ky,tm,tl)
  u0=exp(2*pi*1i.*(kx.*xx(tm)+ky.*yy(tm))); % illumination
  u2=ft(u0.*tm); % lcos illumination
  u4=ft(u2.*tl); % bfp illumination
  u5=ft(u4); % field simple lens instead of microscope)
end

