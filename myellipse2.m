%myellipse2.m
%
%22.2.23
%version 2.1
%@@@ possible enhancements

clc;
clear variables;

%approximate values of ellipse axes from Galton and Hamilton Dickson
%the important thing is actually the ratio of these two measures
%the values below are in ratio 1.87, conforming with H-D's calculations
my_a= sqrt(7); %major ellipse axis dimension
my_b= sqrt(2); %minor ellipse axis dimension

%basic ellipse is  x2/a2 + y2/b2 = 1, centred on (0,0)
%but general ellipse allows rotation (used here) and translation (not used here)

%set up the figure
figure;

%rotation angle alpha will be used in an affine rotation of the ellipse 
%this rotates the major ellipse axis counterclockwise by angle alpha from the horizontal axis
%@@@add input questions here
onealpha=1; %controls (1=yes, 0=no) if we use just one alpha value or more
nconcentric=4; %number of concentric ellipses generated if onealpha=1
%initialise arrays for saving estimated values
tanalpha=zeros(nconcentric, 1);
phi=zeros(nconcentric,1);
theta=zeros(nconcentric,1);
est_r=zeros(nconcentric,1);
%scaler values must start with 1 and end with the largest
%and must have nconcentric values
escaler=[sqrt(2) 0.5 1.5 2]; %the scaler values for the concentric ellipses
if onealpha==1
    %just one alpha value
    myalpha=atan(0.5); %approximate value from Galton and Hamilton Dickson
    %myalpha=pi/4; %rotated midway between axes
    wantaxis=1;
    wanttang=1;
    if nconcentric==1
        mytext={'Galton''s stature ellipse','x=Mid-parent height deviations','y=Offspring height deviations'};
    else
        mytext={'Galton''s stature ellipses','x=Mid-parent height deviations','y=Offspring height deviations'};
    end
    for j=1:nconcentric
        disp(['Ellipse number ' num2str(j)]);
        ellipse_scaler=escaler(j);
        if j==1
            disp('This ellipse was drawn based on Hamilton Dickson Equation 7');
            %draw using H-D equation 7 and given values
            %given values
            HDa=1.22;
            HDb=1.5;
            HDtantheta=2/3;
            %first we have to compute the values for the general ellipse equations
            %this is a check on H-D's working essentially
            %use given values and H-D equations to compute general ellipse coeffs A, B, C etc
            %HD equation 7 implies these values for general ellipse
            A=HDa*HDa;
            B=-2*HDa*HDa*HDtantheta;
            C=HDb*HDb + HDa*HDa*HDtantheta*HDtantheta; %same as HD's liitle c
            D=0; %as elipse centre is (0,0)
            E=0; %as elipse centre is (0,0)
            F=-HDa*HDa*HDb*HDb; %as elipse centre is (0,0)
            %then use A, B etc to calc general ellipse equations values (a,b,alpha) from wiki
            %ellipse page; alpha is shown as theta on wiki page
            %basically using the extended formula for a quadratic as applied to
            %a quadric with terms in y2 xy and x2
            numterm1 = 2*(B*B -4*A*C)*F;
            numterm2a = (A + C + sqrt((A-C)*(A-C) +B*B));
            numterm2b = (A + C - sqrt((A-C)*(A-C) +B*B));
            denomterm=(B*B - 4*A*C);
            calca = -sqrt(numterm1*numterm2a)/denomterm;
            calcb = -sqrt(numterm1*numterm2b)/denomterm;
            calctanalpha=(C - A -sqrt((A-C)*(A-C) +B*B))/B ;
            %now draw the ellipse using these calculated values
            myfig={'-m'}; %more properties to add in cell array
            [x,y,tanalpha(j),phi(j),theta(j)]=drawellipse(ellipse_scaler*calca,ellipse_scaler*calcb,atan(calctanalpha),mytext,myfig,wantaxis,wanttang);
        else 
            disp('This ellipse was drawn using the approximate values');
            disp('given by H-D for major:minor axis ratios and ellipse rotation angle');
            %draw using a standard ellipse equation using 
            %approx values computed and given by H-D in appendix
            myfig={'-b'}; %more properties to add in cell array
            [x,y,tanalpha(j),phi(j),theta(j)]=drawellipse(ellipse_scaler*my_a,ellipse_scaler*my_b,myalpha,mytext,myfig,wantaxis,wanttang);
        end
        %now derive some key values from the drawn ellipses
        [max_y, max_yi]=max(y); %get the maximum y value and the index of this value
        [max_x, max_xi]=max(x);  %get the maximum x value and the index of this value
        ellip_rad=sqrt(y.*y +x.*x); %computes the ellipse radii
        [majorax,maj_index] = max(ellip_rad); %get the major axis (half) length and index
        [minorax,min_index] = min(ellip_rad); %get the minor axis (half) length and index
        est_r(j)=sqrt(phi(j)*theta(j)); %this is one formula for the correlation coeff
        %display the results for each concentric circle -- should be the
        %same values
        disp('Next line has ellipse scaler value, plus computed values of Phi, Theta, and corr');
        disp([ellipse_scaler,phi(j),theta(j),est_r(j)]); 
        %major:minor axis ratios
        if j==1
            %now report the calculated values based on the HD equation 7
            %as the values are rescaled we just report the ratio here
            calc_majminratio=calca/calcb;
            disp(['Major: minor axis ratio used to draw this ellipse = ' num2str(calc_majminratio)]);
        else
            %now report the approximate values given by HD
            %which we used with the general ellipse equation ellipses
            %as the values are rescaled we just report the ratio here
            set_majminratio=my_a/my_b;
            disp(['Major: minor axis ratio used to draw this ellipse = ' num2str(set_majminratio)]);
        end
        if j==1
            set_alpha_degrees=(180/pi)*atan(calctanalpha);
        else
            set_alpha_degrees=(180/pi)*myalpha;
        end
        %now use the end of the major axis point of contact with ellipse to compute alpha angle
        est_alpha_degrees=(180/pi)*atan(y(maj_index)/x(maj_index));
        %use the returned value based on the general ellipse equation
        %drawing routine
        est_alpha_degrees2=(180/pi)*atan(tanalpha(j));
        disp(['Value of ellipse rotation used to draw the ellipse = ' num2str(set_alpha_degrees)]);
        disp(['Computed value of ellipse rotation = ' num2str(est_alpha_degrees)]);
        disp(['Another computed value of ellipse rotation = ' num2str(est_alpha_degrees2)]);
        disp(' ');
        %interesting observation
        %the angle of the line from the origin through the points
        %of intersection of the horiz and vertical tangents
        %is at an angle of atan(1/sqrt(2)) to the x-axis
        %1/sqrt(2) is the value of e/b in Hamilton Dickson's analysis
        yxscaling=y(max_yi)/x(max_xi); %equals 1/sqrt(2)=sqrt(2)/2
    end
elseif onealpha==0
    %range of alpha values
    endalpha=pi/2;
    wantaxis=0;
    mytext={'Rotated ellipses centred on [0,0]','X','Y'};
    myfig={'*-b'}; %more properties to add in cell array
    for myalpha = 0:pi/12:endalpha %pi:-pi/16:pi/2
        drawellipse(my_a,my_b,myalpha,mytext,myfig,wantaxis,[]);
        if myalpha~=endalpha
            hold on;
            disp('Hit key for next ellipse');
            pause
        end
    end

end

function [xt, yt, tanalpha, vgradient, hgradient] = drawellipse(a,b,alpha,figtext,plotprops,doaxis,dotang)

        %alpha in this function is the counter-clockwise angle of the major axis 
        %from the x-axis
        %vary t the angle parameter used to generate the ellipse
        t= 0:0.01:2*pi; %range of eccentricity parameter taken from wikipedia, where it is denoted t
        
        %test for degeneracy of the ellipse
        A = a*a*sin(alpha)*sin(alpha) + b*b*cos(alpha)*cos(alpha);
        B = 2*(b*b - a*a)*sin(alpha)*cos(alpha); %this is non zero when ellipse is rotated
        C = a*a*cos(alpha)*cos(alpha) + b*b*sin(alpha)*sin(alpha);
        D = 0; %for (0,0) centred ellipse
        E = 0; %for (0,0) centred ellipse
        F = -a*a*b*b; %for (0,0) centred ellipse
        mydet = (A*C - B*B/4)*F + B*E*D/4 -C*D*D/4 -A*E*E/4;
        if C*mydet>=0
            disp('This was a degenerate case');
            disp(a); disp(b); disp(alpha);
            return
        end
        if alpha==0
            tanalpha=0;
        elseif alpha==pi/2
            tanalpha=Inf;
        else    
            tanalpha=(C - A -sqrt((A-C)*(A-C) +B*B))/B ;%from wiki ellipse page
        end
        x = a.*cos(t);
        y = b.*sin(t);
        %xt = cos(alpha).*x + sin(alpha).*y; %used with decreasing angle
        xt = cos(alpha).*x - sin(alpha).*y;
        %yt = -sin(alpha).*x + cos(alpha).*y; %used with decreasing angle
        yt = sin(alpha).*x + cos(alpha).*y;
        %draw the ellipse
        plot(xt,yt,char(plotprops(1)),'LineWidth',2.0);
        %plot(x,y,'*-k');
        if doaxis==1
            [max_yt,ymaxi]=max(yt);
            [max_xt,xmaxi]=max(xt);
            plotxmax =ceil(max_xt); %+0.5;
            plotx4ymax=xt(ymaxi);
            plotymax=ceil(max_yt); %+0.5;
            ploty4xmax=yt(xmaxi);
            axmax=max(plotxmax,plotymax);
            plotxmax=axmax;
            plotymax=axmax;
            axis([-plotxmax plotxmax -plotymax plotymax]);
            hold on;
            plot([-plotxmax,plotxmax],[0,0],'-k','LineWidth',1.5);
            hold on;
            plot([0,0],[-plotymax,plotymax],'-k','LineWidth',1.5);
        end    
        if dotang==1
            hold on;
            %plot([0,plotxmax],[max_yt,max_yt],'-r','LineWidth',1.5); %horizontal tangent, plotted past contact
            plot([0,plotx4ymax],[max_yt,max_yt],'-r','LineWidth',1.5); %horizontal tangent, plotted to contact
            hold on;
            plot([0,plotx4ymax],[0,max_yt],'-r','LineWidth',1.5); %line from centre to h tangent intersection
            hold on;
            hgradient=plotx4ymax/max_yt;
            %plot([max_xt,max_xt],[0,plotymax],'-g','LineWidth',1.5); %vertical tangent, plotted past contact
            plot([max_xt,max_xt],[0,ploty4xmax],'-g','LineWidth',1.5); %vertical tangent, plotted to contact
            hold on;
            plot([0,max_xt],[0,ploty4xmax],'-g','LineWidth',1.5); %line from centre to v tangent intersection
            vgradient=ploty4xmax/max_xt;
        end
        xlabel(figtext(2),'fontweight','bold','fontsize',14);
        ylabel(figtext(3),'fontweight','bold','fontsize',14);
        title(figtext(1),'fontweight','bold','fontsize',16);

end