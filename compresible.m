%% Matrius de dades (cas 1D)
P=zeros(1,N+1); T=zeros(1,N+1); V=zeros(1,N+1); ro=zeros(1,N+1);
M=zeros(1,N+1); Sgen=zeros(1,N+1); entr=zeros(1,N+1);
%Variable que se activa una vez se produce el salto de M
salt=false;
%Variable que se activa una vez llegamos al punto donde aparecería onda de
%choque
xoc=false;
%% Cas inicial
P(1)=Pin; T(1)=Tin; V(1)=Vin;
ro(1)=P(1)/(R*T(1));
mins=ro(1)*V(1)*S(1);
M(1)=Min; Sgen(1)=0; entr(1)=0;

%% Comienzan iteraciones
for i=1:(N);
    if i>=choq
        xoc=true;
    end
    % Alteración de los valores para que nunca estemos en M=1.
                if M(i)>=0.97 && salt==false
                    Vch=V(i);
                    Pch=P(i);
                    Tch=T(i);
                    roch=ro(i);
                    Mch=M(i);
                    V(i)=V(i)*1.05;
                    P(i)=P(i)*0.99999;
                    T(i)=T(i)*0.99999;
                    ro(i)=P(i)/(R*T(i));
                    mins=ro(i)*V(i)*S(i);
                    M(i)=V(i)/(sqrt((gamma*R*T(i))));
                    salt=true;
                    entra=true;
                end
    %Valores a la salida supuestos
    P(i+1)=P(i);
    T(i+1)=T(i);
    V(i+1)=V(i);
    ro(i+1)=ro(i);
    % Variable que se activa cuando el VC converge
    bien=false;            
while  bien==false
        %Propiedades medias supuestas
        ro_m=(ro(i+1)+ro(i))/2;
        V_m=(V(i+1)+V(i))/2;
        T_m=(T(i+1)+T(i))/2;
        P_m=ro_m*R*T_m;
        Re=ro_m*V_m*D(i)*2/nu(T_m);
        %% T*
        %PARA GASES
        if Re>2000
            Data=Tref(1,d,T_m,P_m,V_m,Tw,D(i),L,R,S(i));
        else
            Data=Tref(2,d,T_m,P_m,V_m,Tw,D(i),L,R,S(i));
        end
        %Propiedades  a T_ref
        Cpi=Data.cp;
        lambdai=Data.lambda;
        nui=Data.nu;
        alphai=Data.a;
        r=Data.r;
        T_r=Data.T_r;
        %Se calcula "f" a partir del Re calculado
                if (Data.Re<2000)
                    f=16/Data.Re;
                elseif (Data.Re>5*10^3 && Data.Re<3*10^4 && er<=0.0001)
                    f=0.079*Data.Re^-0.25;
                elseif (Data.Re>5*10^3 && Data.Re<3*10^4 && er>0.003 && er<0.005)
                    f=0.096*Data.Re^-0.25;
                elseif (Data.Re>3*10^4 && Data.Re<3*10^6 && er<=0.0001)
                    f=0.046*Data.Re^-0.2;
                elseif (Data.Re>3*10^4 && Data.Re<3*10^6 && er>0.003&& er<0.005)
                    f=0.078*Data.Re^-0.2;
                end   
                        %Churchill General Expression
        %         AS=(2.457*log(1/((7/Re)^0.9+0.27*er)))^16;
        %         BS=(37530/Re)^16;
        %         f=2*((8/Re)^12+1/(AS+BS)^(3/2))^(1/12);
        
        %% Cálculo de coeficientes que resuelven ecuaciones N-S
        Al_m=(Al(i+1)+Al(i))/2;
        Av=mins+f*ro_m*abs(V_m)/4*Al_m*cos(theta(i));
        Bv=(S(i)+Al_m/2*sin(theta(i)));
        Cv=(S(i)+Al_m/2*sin(theta(i)))*P(i)+(mins-f*ro_m*abs(V_m)/4*Al_m*cos(theta(i)))*V(i);
        
        At=mins*cp(T_m)+alphai*Al(i)/2;
        Bt=mins/2+r*alphai*Al(i)/(4*Data.cptref);
        Ct=(mins*cp(T_m)-alphai*Al(i)/2)*T(i)+(mins/2-r*alphai*Al(i)/(4*Data.cptref))*V(i)^2+alphai*Tw*Al(i);
        
        A=Av*At*S(i)-Bv*Bt*mins*R;
        B=Cv*At*S(i);
        C=Bv*Ct*mins*R;
        
        dis=B^2-4*A*C;
        %% Únicamente se prosigue si discriminante es positivo
        if dis>0
            %Se obtienen 2 soluciones. 
            %Revisión de las soluciones mediante entropia generada
            
            positiva=1;% para raíz positiva 
            DataFinal_a = checkSgen(dis,A,B,C,Av,Bv,Cv,At,Bt,Ct,positiva,T(i),P(i),mins,Data,dz,Al(i),S(i),Tw,V(i),f,Data.ro);
            Sgena=DataFinal_a.sgen;

            positiva=0;% para raíz negativa
            DataFinal_b = checkSgen(dis,A,B,C,Av,Bv,Cv,At,Bt,Ct,positiva,T(i),P(i),mins,Data,dz,Al(i),S(i),Tw,V(i),f,Data.ro);
            Sgenb=DataFinal_b.sgen;   
            
            if  Sgena>=0 && Sgenb<0 
                V_i=DataFinal_a.v;
                P_i=DataFinal_a.p;
                T_i=DataFinal_a.T;
                ro_i=DataFinal_a.ro;
                % Se revisa convergenica
                    if abs(P_i-P(i+1))<d && abs(T_i-T(i+1))<d && abs(V_i-V(i+1))<d
                        bien=true;
                    end
                 % Calculo de los valores a la salida del VC
                V(i+1)=V_i;
                P(i+1)=P_i;
                T(i+1)=T_i;
                ro(i+1)=ro_i;
                Sgen(i+1)=Sgena;
                entr(i+1)=entr(i)+DataFinal_a.S;
                gamma=(cp(T(i+1))/R)/(cp(T(i+1))/R-1);
                M(i+1)=V(i+1)/(sqrt((gamma*R*T(i+1))));
           
            elseif Sgena<0 && Sgenb>=0 
                V_i=DataFinal_b.v;
                P_i=DataFinal_b.p;
                T_i=DataFinal_b.T;
                ro_i=DataFinal_b.ro;
                    if abs(P_i-P(i+1))<d && abs(T_i-T(i+1))<d && abs(V_i-V(i+1))<d
                        bien=true;
                    end
                V(i+1)=V_i;
                P(i+1)=P_i;
                T(i+1)=T_i;
                ro(i+1)=ro_i;
                Sgen(i+1)=Sgenb;
                entr(i+1)=entr(i)+DataFinal_b.S;
                gamma=(cp(T(i+1))/R)/(cp(T(i+1))/R-1);
                M(i+1)=V(i+1)/(sqrt((gamma*R*T(i+1))));
            else
                % En caso de flujo crítico, se obtienen dos solcuiones
                % positivas
                
                a=abs(DataFinal_a.v-V(i));% Salto 
                b=abs(DataFinal_b.v-V(i));
                if a<b 
                       if i ~=choq
                           % Se escoge el caso sin salto
                            V_i=DataFinal_a.v;
                            P_i=DataFinal_a.p;
                            T_i=DataFinal_a.T;
                            ro_i=DataFinal_a.ro;
                            if abs(P_i-P(i+1))<d && abs(T_i-T(i+1))<d && abs(V_i-V(i+1))<d
                                bien=true;
                            end
                            V(i+1)=V_i;
                            P(i+1)=P_i;
                            T(i+1)=T_i;
                            ro(i+1)=ro_i;
                            Sgen(i+1)=Sgena;
                            entr(i+1)=entr(i)+DataFinal_a.S;
                            gamma=(cp(T(i+1))/R)/(cp(T(i+1))/R-1);
                            M(i+1)=V(i+1)/(sqrt((gamma*R*T(i+1))));
                       else
                            % Se escoge el caso que genera salto (onda de
                            % choque)
                            V_i=DataFinal_b.v;
                            P_i=DataFinal_b.p;
                            T_i=DataFinal_b.T;
                            ro_i=DataFinal_b.ro;
                            if abs(P_i-P(i+1))<d && abs(T_i-T(i+1))<d && abs(V_i-V(i+1))<d
                                bien=true;
                            end
                            V(i+1)=V_i;
                            P(i+1)=P_i;
                            T(i+1)=T_i;
                            ro(i+1)=ro_i;
                            Sgen(i+1)=Sgenb;
                            entr(i+1)=entr(i)+DataFinal_b.S;
                            gamma=(cp(T(i+1))/R)/(cp(T(i+1))/R-1);
                            M(i+1)=V(i+1)/(sqrt((gamma*R*T(i+1))));
                        end
                   else
                         if i ~=choq
                            V_i=DataFinal_b.v;
                            P_i=DataFinal_b.p;
                            T_i=DataFinal_b.T;
                            ro_i=DataFinal_b.ro;
                                            if abs(P_i-P(i+1))<d && abs(T_i-T(i+1))<d && abs(V_i-V(i+1))<d
                                                bien=true;
                                            end
                            V(i+1)=V_i;
                            P(i+1)=P_i;
                            T(i+1)=T_i;
                            ro(i+1)=ro_i;
                            Sgen(i+1)=Sgenb;
                            entr(i+1)=entr(i)+DataFinal_b.S;
                            gamma=(cp(T(i+1))/R)/(cp(T(i+1))/R-1);
                            M(i+1)=V(i+1)/(sqrt((gamma*R*T(i+1))));
                        else
                            V_i=DataFinal_a.v;
                            P_i=DataFinal_a.p;
                            T_i=DataFinal_a.T;
                            ro_i=DataFinal_a.ro;
                            if abs(P_i-P(i+1))<d && abs(T_i-T(i+1))<d && abs(V_i-V(i+1))<d
                                bien=true;
                            end
                            V(i+1)=V_i;
                            P(i+1)=P_i;
                            T(i+1)=T_i;
                            ro(i+1)=ro_i;
                            Sgen(i+1)=Sgena;
                            entr(i+1)=entr(i)+DataFinal_a.S;
                            gamma=(cp(T(i+1))/R)/(cp(T(i+1))/R-1);
                            M(i+1)=V(i+1)/(sqrt((gamma*R*T(i+1))));
                        end
                  end 
            end
        else
            display('Dis negatiu. Error');
           break;     
        end                          
    end
end




