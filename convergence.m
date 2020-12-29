%Test for convergence of the CCS estimates with increasing series length
%Convergence is given by a decreasing rolling standard deviation 

%Inputs:
%Numpoints: number of lengths to test 

%Outputs:
%C: CCS scores for each length
%r: Slope of best fit line for log-linear rolling standard deviation

function [C,r]=convergence(x,y,lag,tau,dim,numpoints)

if lag>=0
    x=x(lag+1:end);
    y=y(1:end-lag);
else
    x=x(1:end+lag);
    y=y(1-lag:end);
end

T=1+(dim-1)*tau;
Lmin=2*T;
Lmax=length(x);

points=genpoints(numpoints,Lmin,Lmax);
C=zeros(numpoints,2);

for i=1:numpoints
    x_t=x(1:points(i));
    y_t=y(1:points(i));
    [ctemp]=CCS(x_t,y_t,0,tau,dim);
    C(i,1)=ctemp(1); C(i,2)=ctemp(2);
end

Cvar=movvar(C,3);

Bxcy=[ones(size(Cvar(:,1))),points]\log10(Cvar(:,1));
Bycx=[ones(size(Cvar(:,2))),points]\log10(Cvar(:,2));
xcyvarfit=Bxcy(1)+Bxcy(2)*points;
ycxvarfit=Bycx(1)+Bycx(2)*points;

r=[Bxcy(2);Bycx(2)];

if nargout==0
    figure
    subplot(1,2,1)
    plot(points,C(:,1))
    hold on
    plot(points,C(:,2))
    axis([0 Lmax 0 1])
    title('Convergence')
    xlabel('Length')
    ylabel('Rho')
    legend('xcy','ycx')
    
    subplot(1,2,2)
    plot(points,log10(Cvar(:,1)))
    hold on 
    plot(points,log10(Cvar(:,2)))
    hold on 
    plot(points,xcyvarfit)
    hold on 
    plot(points,ycxvarfit)
    title('Var Convergence')
    xlabel('Length')
    ylabel('log(Var)')
    xfstr = "xcyVarFit,r="+num2str(r(1));
    yfstr = "ycxVarFit,r="+num2str(r(2));
    legend('xcyVar','ycxVar',xfstr,yfstr)
end
end
