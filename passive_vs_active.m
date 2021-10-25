clc;
clear all;
close all;

%% simulation parameter
repeat=1000;
n_sample=100:100:400;

passive_err_prob=cell(repeat,1);
active_err_prob=cell(repeat,1);
for i=1:repeat
    passive_err_prob{i,1}=zeros(1,length(n_sample));
    active_err_prob{i,1}=zeros(1,length(n_sample));
end

%% model generation
p=300;
theta=0.05;
p_cross=theta;
model='chain';
[edge,Adj]=adj_generation(model,p);

%% theoretical bound
%See Tandon, Tan, and Zhu 2020
N=length(n_sample);
Kp=-log(1-theta*(1-sqrt(4*theta*(1-theta))));
sigma2=theta*sqrt(4*theta*(1-theta))*exp(Kp);
z=sqrt(theta/(1-theta));
tilf=exp(-Kp*n_sample)./sqrt(2*pi*sigma2*n_sample).*[ones(1,N)+(1-3*sigma2)/(8*sigma2)*ones(1,N)./n_sample];
f=tilf/(1-z).*[ones(1,N)-z*(1+z)/(2*(1-z)^2*sigma2)*ones(1,N)./n_sample];
asy_passive_prob=2*f-tilf;
temp_theo=(p-2)*asy_passive_prob;

%% Structure learning
for rep=1:repeat
    for k=1:length(n_sample) 
        %% sample generation
        if(k==1)
            X=samplegeneration(Adj,theta,n_sample(k));
        else
            X_temp=samplegeneration(Adj,theta,n_sample(k)-n_sample(k-1));
            X=[X;X_temp];
        end

        %% Chow-Liu algorithm
        RHO=X.'*X/n_sample(k);
        RHO(logical(eye(p)))=-inf(1,p);
        CL_adjacency=ChowLiu(RHO);
        CL_err=sum(sum((CL_adjacency~=Adj)));
        CL_err=(CL_err>0);
        passive_err_prob{rep,1}(1,k)=CL_err;

        %% Active-LATHE
        [active_adjacency,tempalpha]=active_LATHE(n_sample(k),p,X,Adj,theta);
        active_err=sum(sum(active_adjacency~=Adj));
        active_err=(active_err>0);
        active_err_prob{rep,1}(1,k)=active_err;
        %ite_alpha(rep,k)=tempalpha;
        

    end
end

passive_err_prob_ave=zeros(1,length(n_sample));
active_err_prob_ave=zeros(1,length(n_sample));
for i=1:repeat
    passive_err_prob_ave=passive_err_prob_ave+passive_err_prob{i,1};
    active_err_prob_ave=active_err_prob_ave+active_err_prob{i,1};
end
passive_err_prob_ave=passive_err_prob_ave/repeat;
active_err_prob_ave=active_err_prob_ave/repeat;



