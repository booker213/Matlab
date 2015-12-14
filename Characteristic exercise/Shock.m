%%Function that draws shocks and expansion fans
function Shock
%%Define Constants
alpha=0.5;
Sr=0.35;
zin=0.7;
psimin=(alpha-1)*zin.^2/4;
xia=0;
xib=-(psimin/Sr);
xic=-(4*psimin)/Sr;
xid=-(3*psimin)/Sr;

%%define storage for functions


%%Define zmap functions for up and down domains
function zup = zup(alpha,zin,psiz)
zup=((1-alpha)*zin+sqrt((((1-alpha)*zin).^2)+4*(1-alpha)*psiz))/(2*(1-alpha));
end

function zdown = zdown(alpha,zin,psiz)
zdown= ((1-alpha)*zin-sqrt((((1-alpha)*zin).^2)+4*(1-alpha)*psiz))/(2*(1-alpha));
end


%% Zinv shock
x=linspace(0,xib);
psi1=0;
zup1=zup(alpha,zin,psi1);

plot(x,zup1,'r');
hold on

%% Up expansion fan
x=linspace(xia,xib);
xi=linspace(xia,xib);
psi2= psimin + Sr*xi;
zup2=zup(alpha,zin,psi2);

plot(x,zup2,'r');
hold on;

for phi=0.5:0.02:1
x=linspace(xia,-(psimin)/(Sr*phi.^2));
xi=linspace(0,-(psimin)/(Sr*phi.^2));
psiup= psimin+Sr*(2*phi-1)*xi;
zupchar=zup(alpha,zin,psiup);
plot(x,zupchar,'g');
hold on;
end

x=linspace(xib,xic);
xi=linspace(xib,xic);
psi3= psimin -Sr*(xi)+2*sqrt(-psimin)*sqrt(Sr*xi);
zup3=zup(alpha,zin,psi3);
plot(x,zup3,'r');
hold on

%down expansion fan
x=linspace(xic,xid);
xi=linspace(xic,xid);
psi4= psimin+ Sr*(xic-xi);
zdown1=zdown(alpha,zin,psi4);
plot(x,zdown1,'b');
hold on ;

for phi=0:0.02:0.5
x=linspace(xic,xic + psimin/(Sr*(1-phi).^2));
xi=linspace(xic,xic + psimin/(Sr*(1-phi).^2));
psidown= psimin+Sr*(2*phi-1)*(xi-xic);
zdownchar=zdown(alpha,zin,psidown);
plot(x,zdownchar,'k');
hold on;
end

x=linspace(xid,xia);
xi=linspace(xid,xia);
psi5 = psimin - Sr*(xic-xi)+2*sqrt(-psimin)*sqrt(Sr*(xic-xi));
zdown2=zdown(alpha,zin,psi5);
plot(x,zdown2,'b');
hold on;

end



