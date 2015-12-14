
function [xr,xexact,ur,uexact]=Characteristics(yr)
%%NB call using  [xfinal,xexact,ufinal,uexact]=Characteristics(0.01)

%%Define Inital Values
yq=0;  dy=yr-yq;
uq=1; xq=1;
aq=uq; bq=1; cq=-2*uq^3;
aj=aq; bj=bq; cj=cq;

for i=1:1000
%%dx/dy characteristic
xr= xq + (aj/bj)*dy;

%%du/dx characteristic
ur= uq+ (cj/aj)*(xr-xq);
%%update constants

aj=0.5*(aq+ur);
cj=0.5*(cq-2*ur^3);
bj=1;
end

%%exact solution
xexact= 0.5*((1+4*dy)^0.5+1);
uexact=(1+4*dy)^(-0.5);


