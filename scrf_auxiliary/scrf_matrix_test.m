clc;clear;close all

%%% define parameters and functions
% define study regions
xrange=[90, 120];
yrange=[20, 50];

% number of previously generated points
n=49;

% generate the locations of the target point and previously generately points
rng(123);
xt=rand(1)*(xrange(2)-xrange(1))+xrange(1);
yt=rand(1)*(yrange(2)-yrange(1))+yrange(1);
xn=rand(n,1)*(xrange(2)-xrange(1))+xrange(1);
yn=rand(n,1)*(yrange(2)-yrange(1))+yrange(1);

figure('color','w')
plot(xn, yn, 'ok','MarkerSize',10);
hold on
plot(xt, yt, '*r','MarkerSize',10);
hold off

% define function used to calculate distance (dist) and spatial correlation (spcc)
% cal_dist=@(x1,y1,x2,y2) ((x1-x2).^2+(y1-y2).^2).^0.5; % just degree distance, not km distance
cal_spcc=@(dist, clen) exp(-dist./clen);


%%% solve the weight calculation equation (C*W=G)
% C is the correlation matrix among n points
% G is the correlation between the target

% Clen values of all points
cmean=400; cstd=100; % the mean value and std of Clen

runtime=1001;
Wall=nan*zeros(n,runtime);
for rr=1:runtime
    if mod(rr,100)==0
       fprintf('runtime %d\n',rr); 
    end
    
    if rr==1
        clen_n=ones(n,1)*cmean; % uniform Clen
        clen_t=cmean; % uniform Clen
    else
        rng('shuffle');
%         clen_n=normrnd(cmean, cstd, [n, 1]); % random Clen
%         clen_t=normrnd(cmean, cstd, 1); % random Clen
        clen_n=gamrnd(50, 10, [n, 1]); % random Clen
        clen_t=gamrnd(50, 10, 1); % random Clen
    end


    C=nan*zeros(n,n);
    for i=1:n
        for j=1:n
            if i==j
                C(i,j)=1;
            elseif i>j
                distij=lldistkm(xn(i),yn(i),xn(j),yn(j));
                ccij1=cal_spcc(distij,clen_n(i));
                ccij2=cal_spcc(distij,clen_n(j));
                C(i,j)=(ccij1+ccij2)/2;
                C(j,i)=C(i,j);
            end
        end
    end


    G=nan*zeros(n,1);
    for i=1:n
        distit=lldistkm(xn(i),yn(i),xt,yt);
        ccij1=cal_spcc(distit,clen_t);
        ccij2=cal_spcc(distit,clen_n(i));
        G(i)=(ccij1+ccij2)/2;
    end


    W=C\G;
    
    Wall(:,rr)=W;
end


figure('color','w')
y1=nanmedian(Wall(:,2:end),2);
y2=prctile(Wall(:,2:end),25,2);
y3=prctile(Wall(:,2:end),75,2);

x=1:length(y1);
x2 = [x, fliplr(x)];
inBetween = [y2', fliplr(y3')];
fill(x2, inBetween, [0.5, 0.5, 0.5]);
hold on
plot(y1,'or','LineWidth',2)
plot(Wall(:,1),'-k','LineWidth',2)
hold off
