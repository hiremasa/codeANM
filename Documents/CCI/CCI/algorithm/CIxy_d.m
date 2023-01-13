function [H]=CIxy_d(x,y)
Hx=CHx_d(x)+CHx_d(y)-CHxy_d(x,y);
H=Hx;
