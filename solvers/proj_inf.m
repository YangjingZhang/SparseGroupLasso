function [y,rr] =  proj_inf(yinput,ld)
y = yinput;
normy = norm(yinput,Inf);
if normy > ld
    y = yinput/normy*ld;
end
rr = (yinput == y);
end