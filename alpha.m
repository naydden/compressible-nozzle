%% Calculo del coeficiente de transferencia de calor
function alphai = alpha( Re_Tref,Gz_Tref,Pr_Tref,lambdai,D)
            if Re_Tref<2000 && Gz_Tref>10
                        C=1.86;
                        m=1/3;
                        n=1/3;
                        K=(D/L)^(1/3)*(nui/nu(Tw))^(0.14);
            elseif Re_Tref<2000 && Gz_Tref<10
                        C=3.66;
                        m=0;
                        n=0;
                        K=1;
            elseif Re_Tref>2000 
                        C=0.023;
                        m=0.8;
                        n=0.4;
                        K=1;
            end
%Nusselt
Nu=C*Re_Tref^m*Pr_Tref^n*K;
%Coeficeinte de calor por convección
alphai=Nu*lambdai/(D);
end

