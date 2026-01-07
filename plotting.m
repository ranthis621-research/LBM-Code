fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Flat Plate\fp_0AoA.mat";
load(fp)



figure
contourf(u_out(500:1500,500:1750,1),'LineStyle','none','LevelStep',0.001)
c = colorbar;
c.Label.String = "U";
colormap jet
axis equal
axis off
set(gca, 'ydir', 'reverse')
clim([-0.01,0.1])



%%
clc
clear all
close all
fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\couette_0_91.mat";
load(fp)

figure
contourf(u_out(:,:,1),'LineStyle','none')
colorbar
colormap turbo


figure
hold on
for x = 1:100:size(u_out(:,:,1),2)
    plot(u_out(:,x,1))
end




%%
clc
clear all
close all
fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\BL_1.mat";
load(fp)

figure
contourf(u_out(:,:,1),'LineStyle','none')
colorbar
colormap turbo


figure
plot(u_out(:,800,1))


%%
fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\output_2DTEST.mat";
load(fp)


%%
fpgif = '2DOSC.gif';
figure
for i = 1:size(u_mat_frames,3)
    if ~all(u_mat_frames(:,:,i)==0)
        contourf(u_mat_frames(:,:,i),'LineStyle','none','LevelStep',0.001)
        c = colorbar;
        c.Label.String = "U";
        colormap jet
        axis equal
        axis off
        set(gca, 'ydir', 'reverse')
        clim([-0.05,0.05])
    
        drawnow
    
        frame = getframe(gcf);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
    
        if i == 1
            imwrite(A,map,fpgif,'gif','LoopCount',Inf,'DelayTime',0.05)
        else
            imwrite(A,map,fpgif,'gif','WriteMode','append','DelayTime',0.05)
        end
    end

end


%%
figure
contourf(u_mat_frames(:,:,20),'LineStyle','none','LevelStep',0.001)
c = colorbar;
c.Label.String = "U";
colormap jet
axis equal
axis off
set(gca, 'ydir', 'reverse')
clim([-0.1,0.1])

%%
fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Boundary Layer\BL_0_79.mat";
load(fp)

xmax = 2;
x = 350;
y = 1:125;

figure
hold on
yindx = find(abs(u_out(y,x,1)/u_out(125,x,1)-0.99)>0.01,1,"first");
plot(u_out(y(1:end-yindx),x,1)/u_out(125,x,1), y(1:end-yindx)/(125-yindx), 'Marker','.')

fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Boundary Layer\BL_0_91.mat";
load(fp)
yindx = find(abs(u_out(y,x,1)/u_out(125,x,1)-0.99)>0.01,1,"first");
plot(u_out(y(1:end-yindx),x,1)/u_out(125,x,1), y(1:end-yindx)/(125-yindx), 'Marker','.')

fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Boundary Layer\BL_1.mat";
load(fp)
yindx = find(abs(u_out(y,x,1)/u_out(125,x,1)-0.99)>0.01,1,"first");
plot(u_out(y(1:end-yindx),x,1)/u_out(125,x,1), y(1:end-yindx)/(125-yindx), 'Marker','.')

fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Boundary Layer\BL_2.mat";
load(fp)
yindx = find(abs(u_out(y,x,1)/u_out(125,x,1)-0.99)>0.01,1,"first");
plot(u_out(y(1:end-yindx),x,1)/u_out(125,x,1), y(1:end-yindx)/(125-yindx), 'Marker','.')

fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Boundary Layer\BL_5.mat";
load(fp)
yindx = find(abs(u_out(y,x,1)/u_out(125,x,1)-0.99)>0.01,1,"first");
plot(u_out(y(1:end-yindx),x,1)/u_out(125,x,1), y(1:end-yindx)/(125-yindx), 'Marker','.')


BL = readmatrix("results\Boundary Layer\blasius_BL_profile.xlsx");
yindx = find(abs(BL(:,2)-0.99)>0.01,1,"last");
plot(BL(:,3),BL(:,1)/5)

legend({"$\tau=0.79$","$\tau=0.91$","$\tau=1$","$\tau=2$","$\tau=5$","Blasius"},'Interpreter','latex')

grid on
xlabel("$U/U_\infty$",'Interpreter','latex')
ylabel("$y$",'Interpreter','latex')


%%
fp = "C:\Users\ranth\OneDrive\Desktop\Projects\LatticeBoltzmannCode\results\Boundary Layer\BL_0_79_larger_domain.mat";
load(fp)


figure
hold on

y = 1:375;
leg_list = {};
n = 1;
for x = ceil(linspace(100,1150,10))

    yindx = find(abs(u_out(y,x,1)/u_out(375,x,1)-0.99)>0.01,1,"first");
    plot(u_out(y(1:end-yindx),x,1)/u_out(375,x,1), y(1:end-yindx)/(375-yindx), 'Marker','.')
    leg_list{n} = "$" + string(x) + "\Delta x$";
    n = n + 1;
end

BL = readmatrix("results\Boundary Layer\blasius_BL_profile.xlsx");
yindx = find(abs(BL(:,2)-0.99)>0.01,1,"last");
plot(BL(:,3),BL(:,1)/5)
leg_list{n} = "Blasius";

legend(leg_list,'Interpreter','latex')

grid on
xlabel("$U/U_\infty$",'Interpreter','latex')
ylabel("$y/\delta$",'Interpreter','latex')


%%
gamma = 1.4;
M = 4;
Merr = 0.05;
dM = Merr/10;
p = 0.1;
perr = 0.006; % psi
dp = perr/10;
p0 = 13.4;
p0err = 0.1; % psi
dp0 = p0err/10;

pi_isentropic = @(M)((1+0.5*(gamma-1)*M.^2).^(-gamma/(gamma-1)));
Cp = @(p,p0,M)((p-pi_isentropic(M)*p0)/(0.5*gamma*pi_isentropic(M)*p0*M^2));


dCp_dp = (Cp(p+dp,p0,M)-Cp(p-dp,p0,M))/(2*dp);
dCp_dp0 = (Cp(p,p0+dp0,M)-Cp(p,p0-dp0,M))/(2*dp0);
dCp_dM = (Cp(p,p0,M+dM)-Cp(p,p0,M-dM))/(2*dM);

Cpcalc = Cp(p,p0,M);
Cperr = sqrt((dCp_dp*perr).^2+(dCp_dp0*p0err).^2+(dCp_dM*Merr).^2);

fprintf("\n\n%f +/- %f \n\n", Cpcalc, Cperr)

%%
a = gcf;

b = a.Children(2).Children;

hold on
for i = 1:length(b)
    xtemp = b(i).XData;
    ytemp = b(i).YData;
    ctemp = b(i).Color;

    p = pi_isentropic(M)*p0*(0.5*gamma*M^2.*ytemp+1);

    dCp_dp = (Cp(p+dp,p0,M)-Cp(p-dp,p0,M))/(2*dp);
    dCp_dp0 = (Cp(p,p0+dp0,M)-Cp(p,p0-dp0,M))/(2*dp0);
    dCp_dM = (Cp(p,p0,M+dM)-Cp(p,p0,M-dM))/(2*dM);

    Cperr = sqrt((dCp_dp*perr).^2+(dCp_dp0*p0err).^2+(dCp_dM*Merr).^2);

    errorbar(xtemp,ytemp,Cperr,Cperr,'vertical','Marker','none','CapSize',8,'LineWidth',1,'Color',ctemp,'LineStyle','none','HandleVisibility','off')
end



%%

a = gcf;