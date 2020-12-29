%Computes number of false nearest neightbors(fnn) across a range of dimensions
%and delays. The appropriate set of parameters should balance relatively
%few dimensions with a low number of FNNs. Simultaneously sweeping over
%both dimension and delay produces a more reliable parameter estimate than
%independently finding each parameter. 

%Krakovská, A., Mezeiová, K., & Budá?ová, H. (2015). Use of false nearest neighbours for selecting variables and embedding parameters for state space reconstruction. Journal of Complex Systems, 2015.

%Input:
%x: T x 1 vector of timepoints 
%maxtau: maximum delay to test 
%maxdim: maximum dimension to test 
%n: number of nearest neighbors to compare. n should be chosen to maximize
%the discrimination between different parameter combinations. Recommended
%range:1-10 

%Output:
%fnn: % FNN as a function of dim and tau
function [fnn]=fnnTDsweep(x,maxtau,maxdim,n)

fnn=zeros(maxdim,maxtau);
for d=1:maxdim
    for tau=1:maxtau
        tic
        dim=d;
        L=length(x);
        T=1+(dim)*tau;
        X=zeros((L-T+1),dim);
 
        
        for t=1:(L-T+1)
            X(t,:)=x((T+t-1):-tau:(T+t-1-(dim-1)*tau));
        end
       
        dX=pdist2(X,X);
        dX=dX+eye(length(dX))*max(max(dX));
        [~,Xs]=sort(dX,2);
        nnX=Xs(:,1:n);
       
        dim=d+1;
        L=length(x);
        T=1+(dim-1)*tau;
        X1=zeros((L-T+1),dim);
        
        for t=1:(L-T+1)
            X1(t,:)=x((T+t-1):-tau:(T+t-1-(dim-1)*tau));
            
        end
        
        dX1=pdist2(X1,X1);
        dX1=dX1+eye(length(dX1))*max(max(dX1));
        [~,Xs]=sort(dX1,2);
        nnX1=Xs(:,1:n);
 
        tnnX=0;
        for i=1:length(dX)
            tnnX=tnnX+numel(intersect(nnX(i,:),nnX1(i,:)));
           
        end
        fnn(d,tau)=1-tnnX/numel(nnX);
    
        toc
    end
end

figure
surf(fnn)
ylabel('Dimension')
xlabel('Tau')
zlabel('% fnn')


end