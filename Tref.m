%% Proceso iterativo para el cálculo de Tref
function Datos = Tref(caso,d,T_m,P_m,V_m,Tw,D,L,R,S)
T_reff=T_m;
T_rn=T_m+1;
listo=false;
while listo==false;
        cp=(1022-0.1626*T_reff+(3.5025*10^-4)*T_reff^2);
        nu=(2.53928*10^-5*sqrt(T_reff/273.15))/(1+(122/T_reff));
        lambda=(3.807*10^-3+7.4*10^-5*T_reff);
        Pr=(cp*nu/lambda);
        %Según turubulento(1) o laminar (2)
        if caso==1
            r=Pr^(1/3);
        else
            r=Pr^(1/2);
        end
        cptref=1/(T_rn-T_m)*(1022*(T_rn-T_m)-0.1626/2*(T_rn^2-T_m^2)+((3.5025*10^-4)/3)*(T_rn^3-T_m^3));
        T_rn=T_m+r*(V_m)^2/(2*cptref);
        T_ref=T_m+0.5*(Tw-T_m)+0.22*(T_rn-T_m);
        
        % Revisamos convergencia
        if abs(T_ref-T_reff)<d
            listo=true;
            T_reff=T_ref;
            T_r=T_rn;
        else
            T_reff=T_ref;
            T_rn=T_rn;
        end
end
ro_Tref=P_m/(R*T_reff);
Ref=ro_Tref*V_m*D/nu;
Gz=pi*D/(4*L)*Ref*Pr;
a=alpha(Ref,Gz,Pr,lambda,D);
% Tupla que devuelve la funcion
Datos = struct('cp',cp,'nu',nu,'lambda',lambda, 'Pr', Pr, 'r', r, 'T_ref',T_reff, 'T_r', T_r, 'Re', Ref, 'Gz', Gz, 'a',a, 'ro', ro_Tref,'cptref',cptref);
end