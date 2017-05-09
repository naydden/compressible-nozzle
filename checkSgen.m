function DataFinal = checkSgen( dis,A,B,C,Av,Bv,Cv,At,Bt,Ct,positiva,Tc,pc,min,Data,dz,Al,S,Tw,V,f,rop)
R=287;   
cp=@(T) 1022-0.1626*T+3.5025*10^-4*T^2;
    if positiva
        v=(B+sqrt(dis))/(2*A);
    else
        v=(B-sqrt(dis))/(2*A);
    end
p=(Cv-Av*v)/Bv;
T=(Ct-Bt*v*v)/At;
ro=p/(R*T);

%Check with conservation equations
% Ec de la masa
ecm=min-ro*v*S;
% Ec del momentum            
eccm=min*(v-V)-pc*S+p*S+f*(ro+rop)/2*abs((v+V)/2)*((v+V)/2)/2*Al;
% Ec de l'energia
ece=min*cp((T+Tc)/2)*(T-Tc)+min*(v^2/2-V^2/2)-Data.a*(Tw-((T+Tc)/2+Data.r*(((v+V)/2)^2)/(2*Data.cptref)))*Al;

%Calculamos Entropia generada           
Ds=1022*log(T/Tc)-0.1626*(T-Tc)+((3.5025*10^-4)/2)*(T^2-Tc^2)-R*log(p/pc);
Sgen=min*Ds/(S*dz)-Data.a*Al*(Tw-Data.T_r)/(S*dz*Tw);
DataFinal = struct('sgen',Sgen,'Ds',Ds,'v',v, 'p', p, 'T', T, 'ro',ro, 'S', Ds); 
end

