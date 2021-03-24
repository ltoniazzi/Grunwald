function [ out ] = pseudosim( input_args )

v=0; D=0.1; pp=0.4; alp=1.9;
pseudo=@(s) 0.1*s.^alp+50*(exp(-pp*s)-1+pp*s)-v*s+D*s.^2; %+10*(exp(-2*pp*s)-1+2*pp*s)+70*(exp(-(1.5)*pp*s)-1+(1.5)*pp*s)
nx=2000;
dt=0.0005;
TFin=1;
BC=4;


%% setup
x=linspace(0,1,nx)';
dx=x(2)-x(1);
M=PseudoOpmatrixwithBC(BC,pseudo,nx,dx,1);



%sim-loop
u0=@(x) exp(-(x-.3).^2*5000)/sqrt(pi/5000);
u=u0(x);
figure(1)
h=plot(x,u);
pause(dt)
for t=dt:dt:TFin
    u=(eye(nx)-dt*M)\u;
    h.YData=u;
    %plot(x,u)
    %loglog(x,u,'-*')
    title(['t = ',num2str(t,'%.4f')])
    ylim([0.0001,5])
    xlim([0,1])
    grid on
    pause(dt)
end

