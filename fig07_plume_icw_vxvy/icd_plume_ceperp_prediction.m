% Simple model of FPC Signature of Ion Cyclotron Damping
%  Greg Howes, 21 JUL 2022

% Modified for lm_fpc
% by Rui Huang, 19 FEB 2025
% iCD Cex Cey signature prediction

close all
clear
clc

save_figure = true;

'======================================================================='
% Set Figure Size
scrsz = get(0, 'ScreenSize');

% data=dlmread('/Users/ghowes/proj/mmsic/icd_case/icd1a_iCD.dat');
data = dlmread('icd1a_iCD.dat');

% icd1a_iCD.dat is generated using the .sm script from plume raw output
% columns in icd1a_iCD.dat are:

% 1 kpar: plume raw output kpar
% 2 gw: normalized total damping rate |\gamma|/|\omega|
% 3 exmag: magnitude of the raw output Ex
% 4 phiex: phase of Ex
% 5 eymag: magnitude of the raw output Ey
% 6 phiey: phase of Ey
% 7 ux1mag: magnitude of the re-normalized (from vA to vTi) Ux
% 8 phiux1: phase of Ux
% 9 uy1mag: magnitude of the re-normalized (from vA to vTi) Uy
% 10 phiuy1 - phiey: phase difference between Uy and Ey

kzval = [0.4 0.5 0.525 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0]
% ikzall = [322 341 345 357 370 382 392 401 417 430 442 452 462]
ikzall = [345] % We only want the 345th row with kpar = 0.525
nkz = size(ikzall, 2) % We only want one kpar value

for ik=1:nkz
    ikz=ikzall(ik)
    
    %theta_loop=(0:15.:180.)
    theta_loop=0.
    ntheta=size(theta_loop,2)
    for itheta=1:ntheta
    
        if (0==1) 
            fac=1.;
            ep1=1.*fac;
            ep2=1.*fac;
            up1=1.*fac;
            up2=1.*fac;
            phi=(0.)*pi;
            d1=(0.49)*pi;
            d2=(0.49)*pi;
        else
            %fac=1.0;
            fac=0.5;
            ep1=data(ikz,3)*fac;
            ep2=data(ikz,5)*fac;
            up1=data(ikz,7)*fac;
            up2=data(ikz,9)*fac;
            phi=data(ikz,6)*pi;
            d1=data(ikz,8)*pi;
            d2=data(ikz,10)*pi;
            kz=data(ikz,1);
            gw=data(ikz,2);
        end
        
        %Set phase shift scan
        alpha=(-pi:pi/20:pi);
        nn=size(alpha,2);
        
        %for in=1:nn
        %d2=alpha(in)
        %d1=d2
        % Species Parameters
        q=1.;
        
        %Create Velocity Space Grid (v_tp units!)
        vmag=4.;
        nvv=64.;
        %Compute vmin, vmax, and dv
        dv=2.*vmag/nvv;
        vmin=-1.*vmag+dv/2.;
        vmax=vmag-dv/2.;
        
        vx=(vmin:dv:vmax)';
        vy=(vmin:dv:vmax)';
        [VX VY]=meshgrid(vx,vy);
        nvx=size(vx,1);
        nvy=size(vy,1);
        
        
        %Time/Phase Variable
        dphi=pi/64;
        omt=(0:dphi:2.*pi-dphi)';
        % omt = (0:dphi:4.*pi)';
        % omt=(0:dphi:8.*pi-dphi)';
        nt=size(omt,1);
        k_times_r = 0.01 * 0.1 + 0.525 * 0.1;
        
        %Create fluid velocity and electric field time series
        b0=1.0;
        u1=zeros(nt);
        u2=zeros(nt);
        e1=zeros(nt);
        e2=zeros(nt);
        % e1=ep1*cos(omt);
        % e2=ep2*cos(omt-phi);
        % u1=up1*cos(omt-d1);
        % u2=up2*cos(omt-phi-d2);
        e1=ep1*cos(k_times_r - omt);
        e2=ep2*cos((k_times_r - omt)+phi);
        u1=up1*cos((k_times_r - omt)+d1);
        u2=up2*cos((k_times_r - omt)+phi+d2);
        
        
        % Rotate by angle theta========================================
        if (1==1)
          theta=(theta_loop(itheta))*pi/180.;
          fprintf('Rotating Wave by angle theta= %7.2f deg or %7.2f rad\n',theta*180/pi,theta);
          e1r=e1; e2r=e2; u1r=u1; u2r=u2;
          e1r=e1*cos(theta) + e2*sin(theta);
          e2r=-e1*sin(theta) + e2*cos(theta);
          u1r=u1*cos(theta) + u2*sin(theta);
          u2r=-u1*sin(theta) + u2*cos(theta);
          e1=e1r;
          e2=e2r;
          u1=u1r;
          u2=u2r;
        end
        % End Rotate by angle theta========================================
        
        % Create variable for ion distribution function, derivatives, and correlations
        f=zeros(nvx,nvy,nt);
        dfdy=zeros(nvx,nvy);
        dfdx=zeros(nvx,nvy);
        cey=zeros(nvx,nvy,nt);
        cex=zeros(nvx,nvy,nt);
        cey_avg=zeros(nvx,nvy);
        cex_avg=zeros(nvx,nvy);
        
        % Time/phase loop========================================================
        for it=1:nt
            %fprintf('Time Loop it= %3i\n',it)
          
            %Create distribution at current time/phase omt: centered at (u1,u2)
            for j=1:nvy
                for i=1:nvx  
                  f(i,j,it)=exp( -((vx(i)-u1(it)).^2. + (vy(j)-u2(it)).^2.));
                end
            end
            %======================================================================
            %Compute v_x and v_y derivatives of f
            dfdy=0;
            dfdx=0;
            for i=1:nvx
                for j=2:nvy-1
                    dfdy(i,j)=(f(i,j+1,it)-f(i,j-1,it))/(vy(j+1)-vy(j-1));
                end
                %End points
                dfdy(i,1)=(f(i,2,it)-f(i,1,it))/(vy(2)-vy(1));
                dfdy(i,nvy)=(f(i,nvy,it)-f(i,nvy-1,it))/(vy(nvy)-vy(nvy-1));
            end
            for j=1:nvy
                for i=2:nvx-1
                    dfdx(i,j)=(f(i+1,j,it)-f(i-1,j,it))/(vx(i+1)-vx(i-1));
                end
                %End points
                dfdx(1,j)=(f(2,j,it)-f(1,j,it))/(vx(2)-vx(1));
                dfdx(nvx,j)=(f(nvx,j,it)-f(nvx-1,j,it))/(vx(nvx)-vx(nvx-1));
            end
          
            %======================================================================
            %Compute FPC 
            for i=1:nvx
                for j=1:nvy
                    if (1==0) % Total v^2
                        cey(i,j,it)=-q*(vx(i)^2+vy(j)^2)/2.*e2(it)*dfdy(i,j);
                        cex(i,j,it)=-q*(vx(i)^2+vy(j)^2)/2.*e1(it)*dfdx(i,j);
                    else % Just component vy^2 or vx^2
                        cey(i,j,it)=-q*vy(j)^2/2.*e2(it)*dfdy(i,j);
                        cex(i,j,it)=-q*vx(i)^2/2.*e1(it)*dfdx(i,j);
                    end
                end
            end
        
          
        end % End time/phase loop=================================================
        
        %Compute the time average
        cex_avg=sum(cex,3)/nt;
        cey_avg=sum(cey,3)/nt;
        ceperp=cex_avg+cey_avg;
        
        %Compute the net energization (from average)
        %NOTE: It is assumed that dvx-dvy here (isotropic bin sizes)
        exsum=sum(cex_avg,'all')*dv*dv;
        eysum=sum(cey_avg,'all')*dv*dv;
        etot=exsum+eysum;
          fprintf('exsum= %8.4f  eysum= %8.4f   etot= %8.4f \n',exsum,eysum,etot);
        
        
        % OUTPUT PLOTS
        %================================================================================
        % PLOTS FOR NAtComm Paper
        %================================================================================
        if (1==0) %===============================================================
          h12=figure('Position',[1 scrsz(4) 2.2*scrsz(3)/3 scrsz(3)/3]);
          ampfac=1E+3;
          
          setauto=1;
          amp1=2.50E-3*ampfac;
          amp2=5.0E-3*ampfac;
          
            t = tiledlayout(1,2);
        %  t.Padding = 'none';
        %  t.Padding = 'compact';
        %  t.TileSpacing = 'none';
        %  t.TileSpacing = 'compact';
        
        %----------------------------------------
        % Tile 1: cey averaged
        ax1=nexttile;
          tmpf(:,:)=cex_avg(:,:);
          % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
          [C,h] = contourf(VX,VY,transpose(tmpf),50);
          set(h,'edgecolor','none');
        set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
        %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
        colorbar;
        
        xlim([-vmag vmag]);
        ylim([-vmag vmag]);
        cL = caxis;
        if (setauto==1)
        caxis([-max(abs(cL)) max(abs(cL))]);
        else
          caxis([-amp1 amp1]);
        end
        colormap(bluewhitered);
        
        daspect([1 1 1]);
        grid on;
        text(-vmag*1.4,vmag*0.85,'(a)','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        text(-vmag*1,vmag*1.15,'$ C_{E_x}(v_x,v_y)$ (arb. units)','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        xlabel('$v_x/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        ylabel('$v_y/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        
        %OUTPUT PARAMETERS
        %text(-1.5*vmag,vmag*1.11,strcat('$E_1= ',num2str(ep1,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        %text(-0.75*vmag,vmag*1.11,strcat('$u_1= ',num2str(up1,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        %text(-1.5*vmag,-vmag*1.35,strcat('$\phi= ',num2str(phi/pi,'%3.2f'),'\pi$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        text(-0.8*vmag,+vmag*0.85,strcat('$\delta_1= ',num2str(d1/pi,'%3.2f'),'\pi$'),'FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        %text(0.5*vmag,-vmag*1.35,strcat('$\delta_1= ',num2str(d1/pi,'%3.2f'),'\pi$'),'FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        %text(-1.45*vmag,-vmag*1.1,strcat('$j_1E_1= ',num2str(exsum,'%3.2g'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold','Rotation',90)
        %text(-1.45*vmag,-vmag*0.,strcat('$j_2E_2= ',num2str(eysum,'%3.2g'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold','Rotation',90)
        
        text(-0.8*vmag,-vmag*0.87,strcat('$\gamma/\omega= ',num2str(gw,'%3.2f'),'$'),'FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        %text(-0.8*vmag,+vmag*0.9,strcat('$\theta= ',num2str(theta*180./pi,'%3.2f'),'^\circ$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        %text(0.4*vmag,-vmag*1.28,strcat('$\int C_{E_\perp 1}= ',num2str(exsum,'%.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        %----------------------------------------
        % Tile 2: cey averaged
        ax3=nexttile;
          tmpf(:,:)=cey_avg(:,:);
          % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.15,50);set(h,'edgecolor','none')
          [C,h] = contourf(VX,VY,transpose(tmpf),50);set(h,'edgecolor','none')
        set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
        %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
        colorbar
        
        xlim([-vmag vmag]);
        ylim([-vmag vmag]);
        cL = caxis;
        if (setauto==1)
        caxis([-max(abs(cL)) max(abs(cL))]);
        else
          caxis([-amp2 amp2]);
        end
        colormap(bluewhitered);
        
        daspect([1 1 1]);
        grid on
        text(-vmag*1.4,vmag*0.85,'(b)','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        text(-vmag*1,vmag*1.15,'$ C_{E_y}(v_x,v_y)$ (arb. units)','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        xlabel('$v_x/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        ylabel('$v_y/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        
        %OUTPUT PARAMETERS
        %text(-1.5*vmag,vmag*1.11,strcat('$E_2= ',num2str(ep2,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        %text(-0.75*vmag,vmag*1.11,strcat('$u_2= ',num2str(up2,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        text(-0.8*vmag,+vmag*0.85,strcat('$\delta_2= ',num2str(d2/pi,'%3.2f'),'\pi$'),'FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        text(-0.8*vmag,-vmag*0.87,strcat('$k_\parallel \rho_i= ',num2str(kz,'%4.3f'),'$'),'FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        
        %text(-1.7*vmag,-vmag*1.28,strcat('$\int C_{E_\perp 2}= ',num2str(eysum,'%.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        %================================================================================
        pngname = strcat('icd1a_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.png');
        epsname = strcat('icd1a_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.eps');
        %pngname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.png');
        %epsname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.eps');
        
        %pngname = strcat('icd1a_ceperp_theta',num2str(theta_loop(itheta),'%03i'),'.png');
        if save_figure
            print(h12,pngname,'-dpng','-r150');
            print(h12,epsname,'-painters','-depsc','-r150');
        end
        
        %================================================================================
        % OUTPUT SOURCE DATA FOR NatComm Paper
        if (1==0)
          fid=fopen('figS2ab.txt','w');
        
          headfmt='%10s %10s %10s %10s\n';
          fprintf(fid,headfmt,'vperp1','vperp2','ceperp1','ceperp2');
          fmt='%10.4f %10.4f %10.4f %10.4f\n';
          for j=1:nvy
            for i=1:nvx
              fprintf(fid,fmt,VX(i,j),VY(i,j),ampfac*cex_avg(j,i),ampfac*cey_avg(j,i));
            end
          end
          fclose(fid);
        
        end %OUTPUT SOURCE DATA FOR NatComm Paper
        %================================================================================
        
        
          
        end %if ==================================================================
        
        %================================================================================
        
        %================================================================================
        % PLOTS of CEperp total
        %================================================================================
        if (1==0)
          h13=figure('Position',[1 scrsz(4) 4.*scrsz(3)/8 6.*scrsz(4)/8]);
          ampfac=1E+3;
          
          setauto=1;
          amp1=2.50E-3*ampfac;
          amp2=5.0E-3*ampfac;
          
        %    t = tiledlayout(1,2);
        %  t.Padding = 'none';
        %  t.Padding = 'compact';
        %  t.TileSpacing = 'none';
        %  t.TileSpacing = 'compact';
        
        %----------------------------------------
        % Tile 1: cey averaged
        %ax1=nexttile;
          tmpf(:,:)=ceperp(:,:);
          [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
        set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
        %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
        colorbar
        
        xlim([-vmag vmag]);
        ylim([-vmag vmag]);
        cL = caxis;
        if (setauto==1)
        caxis([-max(abs(cL)) max(abs(cL))]);
        else
          caxis([-amp1 amp1]);
        end
        colormap(bluewhitered);
        
        daspect([1 1 1]);
        grid on
        text(-vmag*0.3,vmag*1.11,'$ C_{E_{\perp}}(v_{\perp 1},v_{\perp 2})$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        %text(-vmag*0.3,vmag*1.11,'$ C^{(v^2)}_{E_{\perp}}(v_{\perp 1},v_{\perp 2})$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        xlabel('$v_{\perp 1}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        ylabel('$v_{\perp 2}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        
        %OUTPUT PARAMETERS
        text(-1.25*vmag,vmag*1.11,strcat('$E_1= ',num2str(ep1,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        text(-0.75*vmag,vmag*1.11,strcat('$E_2= ',num2str(ep2,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        text(0.5*vmag,vmag*1.11,strcat('$u_1= ',num2str(up1,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        text(1.*vmag,vmag*1.11,strcat('$u_2= ',num2str(up2,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        text(-0.75*vmag,-vmag*1.25,strcat('$\phi= ',num2str(phi/pi,'%3.2f'),'\pi$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        text(-1.25*vmag,-vmag*1.25,strcat('$k_\parallel \rho_i= ',num2str(kz,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        text(0.25*vmag,-vmag*1.25,strcat('$\delta_1= ',num2str(d1/pi,'%3.2f'),'\pi$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        text(0.85*vmag,-vmag*1.25,strcat('$\delta_2= ',num2str(d2/pi,'%3.2f'),'\pi$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        
        text(-1.25*vmag,-vmag*0.75,strcat('$j_1E_1= ',num2str(exsum,'%3.2g'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold','Rotation',90)
        text(-1.25*vmag,vmag*0.25,strcat('$j_2E_2= ',num2str(eysum,'%3.2g'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold','Rotation',90)
        
        text(-0.8*vmag,-vmag*0.9,strcat('$\gamma/\omega= ',num2str(gw,'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        text(-0.95*vmag,0.9*vmag,'(a)','FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        
        %OUTPUT PARAMETERS
        
        
        
        %text(-1.7*vmag,-vmag*1.28,strcat('$\int C_{E_\perp 2}= ',num2str(eysum,'%.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold')
        
        %================================================================================
        pngname = strcat('icd1a_ceperp',num2str(round(kz*10),'%02i'),'.png');
        epsname = strcat('icd1a_ceperp',num2str(round(kz*10),'%02i'),'.eps');
        %pngname = strcat('icd1a_fullv_ceperp',num2str(round(kz*10),'%02i'),'.png');
        %epsname = strcat('icd1a_fullv_ceperp',num2str(round(kz*10),'%02i'),'.eps');
        
        print(h13,pngname,'-dpng','-r150')
        print(h13,epsname,'-painters','-depsc','-r150')
        
        
        end
        %================================================================================
        % PLOTS of Different ICW wave phases and how they sum
        %================================================================================
        if (1==0)
          h13=figure('Position',[1 scrsz(4) 8.*scrsz(3)/8 6.*scrsz(4)/8]);
          ampfac=1E+2;
          
          setauto=0;
          amp1=20;
        %  amp1=2.50E-3*ampfac;
          amp2=5.0E-3*ampfac;
          
            t = tiledlayout(2,4);
        %  t.Padding = 'none';
          t.Padding = 'compact';
        %  t.TileSpacing = 'none';
          t.TileSpacing = 'compact';
        
        %----------------------------------------
        % Tiles 0 through 7
          for it=0:7
            jt=it*16+8;  %Set index out of 128 over 2 pi
           if (it==0)
             ax1=nexttile;
             elseif (it==1)
             ax2=nexttile;
             elseif (it==2)
             ax3=nexttile;
             elseif (it==3)
             ax4=nexttile;
             elseif (it==4)
             ax5=nexttile;
             elseif (it==5)
             ax6=nexttile;
             elseif (it==6)
             ax7=nexttile;
             elseif (it==7)
             ax8=nexttile;
        end
          tmpf(:,:)=cex(:,:,jt);
          [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
        set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
        
        xlim([-vmag vmag]);
        ylim([-vmag vmag]);
        cL = caxis;
        if (setauto==1)
        caxis([-max(abs(cL)) max(abs(cL))]);
        else
          caxis([-amp1 amp1]);
        end
        colormap(bluewhitered);
        
        
        if (it==3 || it==7)
        c=  colorbar('eastoutside','FontSize',16);
        %  c=  colorbar('eastoutside')
        %c=  colorbar('eastoutside','Position',[0.9505 0.1211 0.0052 0.3663]);
        %set(c,'Position',[0.9805 0.1211 0.0052 0.3663])
           %  c=  colorbar('manual','Position',[0.8566 0.1100 0.0052 0.3885]);
        %  [left, bottom, width, height].
        %  c.Position=[0.5 0.5 0.5 0.5]
          
        end
        
        
        daspect([1 1 1]);
        grid on
        text(-vmag*0.75,vmag*0.8,strcat('$ \omega t/2\pi = ',num2str(omt(jt)/(2.*pi),'%4.3f'),'$'),'Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'Interpreter','latex','FontWeight','bold')
        
        if (it>=4)
          xlabel('$v_{\perp 1}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
          end
        if (it==0 || it==4)
        ylabel('$v_{\perp 2}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',32,'FontWeight','bold')
        end
         end
        
        %================================================================================
        pngname = strcat('icd1a_ceperp1_phases.png');
        epsname = strcat('icd1a_ceperp1_phases.eps');
        %pngname = strcat('icd1a_fullv_ceperp',num2str(round(kz*10),'%02i'),'.png');
        %epsname = strcat('icd1a_fullv_ceperp',num2str(round(kz*10),'%02i'),'.eps');
        
        print(h13,pngname,'-dpng','-r150')
        print(h13,epsname,'-painters','-depsc','-r150')
        
        
        end
        
        
        %================================================================================
        % PLOTS of Ex, Ey, Ux, Uy on separate panels
        %================================================================================
        if (1==0)
            h14=figure('Position',[1 scrsz(4) scrsz(3) 2.*scrsz(3)/8]);

            t2 = tiledlayout(1, 2);
            
            ax21 = nexttile;
            plot(omt(:)/(2.*pi),u1(:),'LineWidth',4.0,'Color','k', 'LineStyle',':', 'DisplayName', "$U_x$");
            hold on;
            plot(omt(:)/(2.*pi),e1(:),'LineWidth',4.0,'Color','k', 'LineStyle','-', 'DisplayName', "$E_x$");
            max_abs_e1_idx = find(abs(e1)==max(abs(e1), [], "all"));
            fprintf("Maximum |Ex| indices are %1.1d, %1.1d \n", max_abs_e1_idx(1), max_abs_e1_idx(2));
            xline(omt(max_abs_e1_idx)/(2.*pi), 'LineWidth',4.0,'Color','k', 'LineStyle','--', 'HandleVisibility','off');

            set(gca,'FontSize',16,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
            legend('Interpreter','latex');
            xlabel('$\Omega t/2\pi$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold')
            xL = xlim;
            yL = ylim;
            line(xL, [0 0],'LineWidth',1.0,'Color','k', 'HandleVisibility','off');  %x-axis
            
            grid on 
            xticks([0. 0.25 0.5 0.75 1]);
            
            ax22 = nexttile;
            plot(omt(:)/(2.*pi),u2(:),'LineWidth',4.0,'Color','k', 'LineStyle',':', 'DisplayName', "$U_y$");
            hold on;
            plot(omt(:)/(2.*pi),e2(:),'LineWidth',4.0,'Color','k', 'LineStyle','-', 'DisplayName', "$E_y$");
            max_abs_e2_idx = find(abs(e2)==max(abs(e2), [], "all"));
            fprintf("Maximum |Ey| indices are %1.1d, %1.1d \n", max_abs_e2_idx(1), max_abs_e2_idx(2));
            xline(omt(max_abs_e2_idx)/(2.*pi), 'LineWidth',4.0,'Color','k', 'LineStyle','--', 'HandleVisibility','off');

            set(gca,'FontSize',16,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
            legend('Interpreter','latex');
            xlabel('$\Omega t/2\pi$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold')
            xL = xlim;
            yL = ylim;
            line(xL, [0 0],'LineWidth',1.0,'Color','k', 'HandleVisibility','off'); %x-axis
            
            grid on 
            xticks([0. 0.25 0.5 0.75 1]);
            

            %================================================================================
            pngname = strcat('icd1a_u1e1u2e2_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.png');
            epsname = strcat('icd1a_u1e1u2e2_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.eps');
            %pngname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.png');
            %epsname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.eps');
            
            %pngname = strcat('icd1a_ceperp_theta',num2str(theta_loop(itheta),'%03i'),'.png');
            if save_figure
                print(h14,pngname,'-dpng','-r150');
                print(h14,epsname,'-painters','-depsc','-r150');
            end
        
        end
        %================================================================================
        % PLOTS of Ex, Ey, Ux, Uy on one figure
        %================================================================================
        if (1==1)
            h16=figure('Position',[1 scrsz(4) scrsz(3) 2.*scrsz(3)/8]);

            % t4 = tiledlayout(1, 1);
            
            % ax21 = nexttile;
            axes("Position", [0.05, 0.28, 0.9, 0.6]);
            plot(omt(:)/(2.*pi),u1(:),'LineWidth',4.0,'Color','k', 'LineStyle',':', 'DisplayName', "$U_x/v_{ti}$");
            hold on;
            plot(omt(:)/(2.*pi),e1(:),'LineWidth',4.0,'Color','k', 'LineStyle','-', 'DisplayName', "$E_x/(v_{ti} B_0 / c)$");
            max_abs_e1_idx = find(abs(e1)==max(abs(e1), [], "all"));
            fprintf("Maximum |Ex| indices are %1.1d, %1.1d \n", max_abs_e1_idx(1), max_abs_e1_idx(2));
            xline(omt(max_abs_e1_idx)/(2.*pi), 'LineWidth',4.0,'Color','k', 'LineStyle','--', 'HandleVisibility','off');

            plot(omt(:)/(2.*pi),u2(:),'LineWidth',4.0,'Color','#c1272d', 'LineStyle',':', 'DisplayName', "$U_y/v_{ti}$");
            hold on;
            plot(omt(:)/(2.*pi),e2(:),'LineWidth',4.0,'Color','#c1272d', 'LineStyle','-', 'DisplayName', "$E_y/(v_{ti} B_0 / c)$");
            max_abs_e2_idx = find(abs(e2)==max(abs(e2), [], "all"));
            fprintf("Maximum |Ey| indices are %1.1d, %1.1d \n", max_abs_e2_idx(1), max_abs_e2_idx(2));
            xline(omt(max_abs_e2_idx)/(2.*pi), 'LineWidth',4.0,'Color','#c1272d', 'LineStyle','--', 'HandleVisibility','off');

            % set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
            % legend('Interpreter','latex');
            % xlabel('$\Omega t/2\pi$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold')
            % xL = xlim;
            % yL = ylim;
            % line(xL, [0 0],'LineWidth',1.0,'Color','k', 'HandleVisibility','off');  %x-axis
            % 
            % grid on 
            % xticks([0. 0.25 0.5 0.75 1]);

            text(0, 0.98, "$(a)$",'FontSize',24,'Interpreter','latex','FontWeight','bold', ...
                'VerticalAlignment', 'bottom');

            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2)
            % legend('Interpreter','latex');
            xlabel('$\Omega t/2\pi$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold')
            xL = xlim;
            yL = ylim;
            line(xL, [0 0],'LineWidth',1.0,'Color','k', 'HandleVisibility','off'); %x-axis
            
            grid on 
            xticks([0. 0.25 0.5 0.75 1]);

            % Add a legend above all subplots
            % hLegend = legend(ax1, 'show');
            hLegend = legend('Interpreter','latex', 'Location','eastoutside', ...
                'Orientation', 'vertical', 'EdgeColor', 'none');
            % set(hLegend, 'Position', [0.5, 0, 1, 1], 'Orientation', 'horizontal');
            % set(hLegend, 'Position', [], 'Orientation', 'horizontal');
            

            %================================================================================
            pngname = strcat('./png/icd1a_u1e1u2e2Single_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.png');
            epsname = strcat('./eps/icd1a_u1e1u2e2Single_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.eps');
            pdfname = strcat('./pdf/icd1a_u1e1u2e2Single_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.pdf');
            %pngname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.png');
            %epsname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.eps');
            
            %pngname = strcat('icd1a_ceperp_theta',num2str(theta_loop(itheta),'%03i'),'.png');
            if save_figure
                print(h16,pngname,'-dpng','-r150');
                print(h16,epsname,'-painters','-depsc','-r150');
                % Match the figure size (convert pixels to inches, or set directly)
                fig_width = 9;   % in inches
                fig_height = 3;   % in inches
                set(h16, 'PaperSize', [fig_width fig_height]);
                set(h16, 'PaperPosition', [0 0 fig_width fig_height]);
                % set(h16, 'PaperPositionMode','Auto');
                print(h16,pdfname,'-dpdf','-vector');
            end
        
        end

        %================================================================================
        % PLOTS of instantaneous Cex and Cey at max |Ex| and |Ey| slices
        %================================================================================
        if (1==1)
            h15=figure('Position',[1 scrsz(4) 6.*scrsz(3)/8 4.*scrsz(3)/8]);
            subpanel_labels = ["$(b)$", "$(c)$", "$(d)$", "$(e)$", "$(f)$", "$(g)$"];

            t3 = tiledlayout(2, 3);
            max_abs_e1_idx = find(abs(e1)==max(abs(e1), [], "all"));
            fprintf("Maximum |Ex| indices are %1.1d, %1.1d \n", max_abs_e1_idx(1), max_abs_e1_idx(2));

            max_abs_e2_idx = find(abs(e2)==max(abs(e2), [], "all"));
            fprintf("Maximum |Ey| indices are %1.1d, %1.1d \n", max_abs_e2_idx(1), max_abs_e2_idx(2));

            %==============================================================
            ax31 = nexttile;
            % tmpf(:,:)=cex(:,:, 2) + cex(:,:, 66);
            tmpf(:, :) = cex(:, :, max_abs_e1_idx(1));
            % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
            [C,h] = contourf(VX,VY,transpose(tmpf),50);
            set(h,'edgecolor','none');
            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
            %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
            colorbar;
            
            xlim([-vmag vmag]);
            ylim([-vmag vmag]);
            cL = caxis;
            caxis([-max(abs(cL)) max(abs(cL))]);
            % if (setauto==1)
            % caxis([-max(abs(cL)) max(abs(cL))]);
            % else
            % caxis([-amp1 amp1]);
            % end
            colormap(bluewhitered);
            
            daspect([1 1 1]);
            grid on;
            text(-0.98*vmag,vmag*1.2,subpanel_labels(1),'Interpreter','latex','FontName','TimesNewRoman','FontSize',38,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-vmag*0,vmag*1,'$C_{E_x}(v_x,v_y)$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-0.9*vmag,vmag*0.87,strcat('$\omega t / (2 \pi)= ',num2str(omt(max_abs_e1_idx(1))/(2*pi),'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold');
            xlabel('$v_{x}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            ylabel('$v_{y}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');

            % Add circle at local bulk flow velocity
            if (1==1)
              hold on;
              ux=u1(max_abs_e1_idx(1)); uy=u2(max_abs_e1_idx(1));
              plot(ux,uy,'LineWidth',2.0,'Color','k','Marker','p','MarkerSize',12);
              alpha=0:2*pi/100:2*pi;
              xp=cos(alpha)+ux;
              yp=sin(alpha)+uy;
              plot(xp,yp,'LineWidth',2.0,'Color','k');
            end

            %==============================================================
            ax32 = nexttile;
            % tmpf(:,:)=cex(:,:, 2) + cex(:,:, 66);
            tmpf(:, :) = cex(:, :, max_abs_e1_idx(2));
            % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
            [C,h] = contourf(VX,VY,transpose(tmpf),50);
            set(h,'edgecolor','none');
            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
            %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
            colorbar;
            
            xlim([-vmag vmag]);
            ylim([-vmag vmag]);
            cL = caxis;
            caxis([-max(abs(cL)) max(abs(cL))]);
            % if (setauto==1)
            % caxis([-max(abs(cL)) max(abs(cL))]);
            % else
            % caxis([-amp1 amp1]);
            % end
            colormap(bluewhitered);
            
            daspect([1 1 1]);
            grid on;
            text(-0.98*vmag,vmag*1.2,subpanel_labels(2),'Interpreter','latex','FontName','TimesNewRoman','FontSize',38,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-vmag*0,vmag*1,'$C_{E_x}(v_x,v_y)$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-0.9*vmag,vmag*0.87,strcat('$\omega t / (2 \pi)= ',num2str(omt(max_abs_e1_idx(2))/(2*pi),'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold');
            xlabel('$v_{x}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            ylabel('$v_{y}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');

            % Add circle at local bulk flow velocity
            if (1==1)
              hold on;
              ux=u1(max_abs_e1_idx(2)); uy=u2(max_abs_e1_idx(2));
              plot(ux,uy,'LineWidth',2.0,'Color','k','Marker','p','MarkerSize',12);
              alpha=0:2*pi/100:2*pi;
              xp=cos(alpha)+ux;
              yp=sin(alpha)+uy;
              plot(xp,yp,'LineWidth',2.0,'Color','k');
            end
            %==============================================================
            ax33 = nexttile;
            % tmpf(:,:)=cex(:,:, 2) + cex(:,:, 66);
            tmpf(:, :) = (cex(:, :, max_abs_e1_idx(1)) + cex(:,:, max_abs_e1_idx(2)))/1.;
            % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
            [C,h] = contourf(VX,VY,transpose(tmpf),50);
            set(h,'edgecolor','none');
            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
            %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
            colorbar;
            
            xlim([-vmag vmag]);
            ylim([-vmag vmag]);
            cL = caxis;
            caxis([-max(abs(cL)) max(abs(cL))]);
            % if (setauto==1)
            % caxis([-max(abs(cL)) max(abs(cL))]);
            % else
            % caxis([-amp1 amp1]);
            % end
            colormap(bluewhitered);
            
            daspect([1 1 1]);
            grid on;
            text(-0.98*vmag,vmag*1.2,subpanel_labels(3),'Interpreter','latex','FontName','TimesNewRoman','FontSize',38,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-vmag*0,vmag*1,'Sum','Interpreter','none','FontName','TimesNewRoman','FontSize',24,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            % text(-0.8*vmag,-vmag*0.87,strcat('$\omega t / 2 \pi= ',num2str(omt(66)/(2*pi),'%3.2f'),'$'),'FontSize',20,'Interpreter','latex','FontWeight','bold');
            xlabel('$v_{x}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            ylabel('$v_{y}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');

            %==============================================================
            ax34 = nexttile;
            % tmpf(:,:)=cey(:,:, 34) + cey(:,:, 98);
            tmpf(:, :) = cey(:, :, max_abs_e2_idx(1));
            % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
            [C,h] = contourf(VX,VY,transpose(tmpf),50);
            set(h,'edgecolor','none');
            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
            %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
            colorbar;
            
            xlim([-vmag vmag]);
            ylim([-vmag vmag]);
            cL = caxis;
            caxis([-max(abs(cL)) max(abs(cL))]);
            % if (setauto==1)
            % caxis([-max(abs(cL)) max(abs(cL))]);
            % else
            % caxis([-amp1 amp1]);
            % end
            colormap(bluewhitered);
            
            daspect([1 1 1]);
            grid on;
            text(-0.98*vmag,vmag*1.2,subpanel_labels(4),'Interpreter','latex','FontName','TimesNewRoman','FontSize',38,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-vmag*0,vmag*1,'$C_{E_y}(v_x,v_y)$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-0.9*vmag,vmag*0.87,strcat('$\omega t / (2 \pi)= ',num2str(omt(max_abs_e2_idx(1))/(2*pi),'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold');
            xlabel('$v_{x}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            ylabel('$v_{y}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            
            % Add circle at local bulk flow velocity
            if (1==1)
              hold on;
              ux=u1(max_abs_e2_idx(1)); uy=u2(max_abs_e2_idx(1));
              plot(ux,uy,'LineWidth',2.0,'Color','k','Marker','p','MarkerSize',12);
              alpha=0:2*pi/100:2*pi;
              xp=cos(alpha)+ux;
              yp=sin(alpha)+uy;
              plot(xp,yp,'LineWidth',2.0,'Color','k');
            end
            %==============================================================
            ax35 = nexttile;
            % tmpf(:,:)=cex(:,:, 2) + cex(:,:, 66);
            tmpf(:, :) = cey(:, :, max_abs_e2_idx(2));
            % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
            [C,h] = contourf(VX,VY,transpose(tmpf),50);
            set(h,'edgecolor','none');
            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
            %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
            colorbar;
            
            xlim([-vmag vmag]);
            ylim([-vmag vmag]);
            cL = caxis;
            caxis([-max(abs(cL)) max(abs(cL))]);
            % if (setauto==1)
            % caxis([-max(abs(cL)) max(abs(cL))]);
            % else
            % caxis([-amp1 amp1]);
            % end
            colormap(bluewhitered);
            
            daspect([1 1 1]);
            grid on;
            text(-0.98*vmag,vmag*1.2,subpanel_labels(5),'Interpreter','latex','FontName','TimesNewRoman','FontSize',38,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-vmag*0,vmag*1, '$C_{E_y}(v_x,v_y)$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-0.9*vmag,vmag*0.87,strcat('$\omega t / (2 \pi)= ',num2str(omt(max_abs_e2_idx(2))/(2*pi),'%3.2f'),'$'),'FontSize',24,'Interpreter','latex','FontWeight','bold');
            xlabel('$v_{x}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            ylabel('$v_{y}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            
            % Add circle at local bulk flow velocity
            if (1==1)
              hold on;
              ux=u1(max_abs_e2_idx(2)); uy=u2(max_abs_e2_idx(2));
              plot(ux,uy,'LineWidth',2.0,'Color','k','Marker','p','MarkerSize',12);
              alpha=0:2*pi/100:2*pi;
              xp=cos(alpha)+ux;
              yp=sin(alpha)+uy;
              plot(xp,yp,'LineWidth',2.0,'Color','k');
            end
            %==============================================================
            ax36 = nexttile;
            % tmpf(:,:)=cex(:,:, 2) + cex(:,:, 66);
            tmpf(:, :) = (cey(:, :, 34) + cey(:,:, 98))/1.;
            % [C,h] = contourf(VX,VY,transpose(tmpf*ampfac)+0.1,50);set(h,'edgecolor','none')
            [C,h] = contourf(VX,VY,transpose(tmpf),50);
            set(h,'edgecolor','none');
            set(gca,'FontSize',24,'FontName','TimesNewRoman','FontWeight','normal','LineWidth',2);
            %colorbar('east','FontSize',16,'FontName','TimesNewRoman')
            colorbar;
            
            xlim([-vmag vmag]);
            ylim([-vmag vmag]);
            cL = caxis;
            caxis([-max(abs(cL)) max(abs(cL))]);
            % if (setauto==1)
            % caxis([-max(abs(cL)) max(abs(cL))]);
            % else
            % caxis([-amp1 amp1]);
            % end
            colormap(bluewhitered);
            
            daspect([1 1 1]);
            grid on;
            text(-0.98*vmag,vmag*1.2,subpanel_labels(6),'Interpreter','latex','FontName','TimesNewRoman','FontSize',38,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            text(-vmag*0,vmag*1,'Sum','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'Interpreter','latex','FontWeight','bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            % text(-0.8*vmag,-vmag*0.87,strcat('$\omega t / 2 \pi= ',num2str(omt(66)/(2*pi),'%3.2f'),'$'),'FontSize',20,'Interpreter','latex','FontWeight','bold');
            xlabel('$v_{x}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            ylabel('$v_{y}/v_{ti}$','Interpreter','latex','FontName','TimesNewRoman','FontSize',24,'FontWeight','bold');
            
            %================================================================================
            pngname = strcat('./png/icd1a_OmtSlices_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.png');
            epsname = strcat('./eps/icd1a_OmtSlices_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.eps');
            pdfname = strcat('./pdf/icd1a_OmtSlices_kz',num2str(round(kz*1000),'%04i'), '_fac', num2str(round(fac*10), "%02i"), '.pdf');
            %pngname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.png');
            %epsname = strcat('icd1a_fullv_kz',num2str(round(kz*10),'%02i'),'.eps');
            
            %pngname = strcat('icd1a_ceperp_theta',num2str(theta_loop(itheta),'%03i'),'.png');
            if save_figure
                print(h15,pngname,'-dpng','-r150');
                print(h15,epsname,'-painters','-depsc','-r150');
                % Match the figure size (convert pixels to inches, or set directly)
                fig_width = 14.4;   % in inches
                fig_height = 9.6;   % in inches
                set(h15, 'PaperSize', [fig_width fig_height]);
                set(h15, 'PaperPosition', [0 0 fig_width fig_height]);
                print(h15,pdfname,'-dpdf','-vector');
            end

        end
    
    
    end % theta loop
    %pause
    %end % Parameter sweep (in)
end % kz sweep (ik)


