function  graficas( DataFINAL,L,N)
figure1 = figure('Color',[1 1 1]);
% legendInfo{ra} = ['Min = ' num2str(Min)]; 
x=linspace(0,L,N+1);
                subplot(3,2,1)
                plot(x,DataFINAL.Pnc,x,DataFINAL.Pc);
                grid on;
                title('Pressió','FontWeight','bold','FontSize',14); 
                xlabel('X [m]','FontWeight','bold'); ylabel('P [Pa]','FontWeight','bold');

                subplot(3,2,2)
                plot(x,DataFINAL.Tnc,x,DataFINAL.Tc);
                hold on;
                grid on;
                title('Temperatura','FontWeight','bold','FontSize',14); 
                xlabel('X [m]','FontWeight','bold'); ylabel('T [K]','FontWeight','bold');

                subplot(3,2,3)
                plot(x,DataFINAL.Vnc,x,DataFINAL.Vc);
                hold on;
                grid on;
                title('Velocitat','FontWeight','bold','FontSize',14); 
                xlabel('X [m]','FontWeight','bold'); ylabel('V [m/s]','FontWeight','bold');

                subplot(3,2,4)
                plot(x,DataFINAL.Mnc,x,DataFINAL.Mc);
                hold on;
                grid on;
                title('Mach','FontWeight','bold','FontSize',14); 
                xlabel('X [m]','FontWeight','bold'); ylabel('Mach','FontWeight','bold');

                subplot(3,2,5)
                semilogy(x,abs(DataFINAL.entrnc),x,abs(DataFINAL.entrc));
                hold on;
                grid on;
                title('Entropia Específica','FontWeight','bold','FontSize',14); 
                xlabel('X [m]','FontWeight','bold'); ylabel('s [J/KgK]','FontWeight','bold');

                subplot(3,2,6)
                plot(x,DataFINAL.Sgennc,x,DataFINAL.Sgenc);
                hold on;
                grid on;
                title('Entropía generada per unitat de Volum','FontWeight','bold','FontSize',14); 
                xlabel('X [m]','FontWeight','bold'); ylabel('W [W/m3K]','FontWeight','bold');
end

