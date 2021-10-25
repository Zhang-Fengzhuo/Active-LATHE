function [active_adjacency,alpha]=active_LATHE(n,p,X,Adj,theta)

%Active-LATHE algorithm to learn tree structure
%Parameters:
%n: number of vector samples
%p: number of nodes
%X: n*p data matrix
%Adj: nominal adjacency matrix of the tree
%theta: the parameter of the homogeneous tree
%
%Output:
%active_adjacency: estimate of the Adjacency matrix
%alpha: the ratio of the number of vector samples used in the Global
%Learning Phase
%
%Fengzhuo Zhang, Oct 2021, NUS


rho_threhold=[0 0.02 0.07 0.16 0.34 0.53 0.76 1];
opt_alpha=[1 0.995 0.985 0.95 0.9 0.85 0.8];

% Phase 1
alpha_i=0.8;    %alpha_{i}
alpha_i1=0.8;   %alpha_{i+1}
while(1)
    alpha_i=alpha_i1;
    num_p1=floor(alpha_i*n);
    X1=X(1:num_p1,:);
    RHO1=X1.'*X1/num_p1;
    RHO1(logical(eye(p)))=-inf(1,p);
    temp_adjacency=ChowLiu(RHO1);
    RHO1(logical(eye(p)))=zeros(1,p);
    rho=sum(sum(RHO1.*temp_adjacency))/2/(p-1);
    indicator=(rho>=rho_threhold);
    indicator=find(indicator);
    indicator=max(indicator);
    indicator=max([1 indicator]);
    alpha_i1=opt_alpha(indicator); %alpha_{i+1}
    if(alpha_i1<=alpha_i)
        break;
    end
end
alpha=alpha_i;

%Phase 2
triple=find_triple(temp_adjacency);
unconfid_triple=find_unconfid_triple(triple,RHO1,rho);

% Learn the triples that share a common edge jointly in unconfident triples
[CP,n_CP]=find_commedge_component(unconfid_triple);
active_adjacency=temp_adjacency;
[n_unconfid_tri,~]=size(unconfid_triple);
if(n_unconfid_tri>0)
    vertices=reshape(unconfid_triple,1,[]);
    vertices=unique(vertices);
    n_vertices=length(vertices);
    num_p2=floor((1-alpha)*n*p/n_vertices);
    X_new=samplegeneration(Adj,theta,num_p1+num_p2-n);
    X_joint=[X;X_new];
    for i=1:n_CP
        index=sort(CP{i,1});
        n_index=length(index);
        X_temp=X_joint(:,index);
        %RHO_temp= X_temp.'*X_temp/(num_p1+num_p2);
        %
        MI_temp=zeros(n_index,n_index);
        for j=1:n_index-1
            for k=j+1:n_index
                data=[X_temp(:,j) X_temp(:,k)];
                MI_temp(j,k)=bi_MI_sample(data);
                MI_temp(k,j)=MI_temp(j,k);
            end
        end
        %}
        subgraph=ChowLiu(MI_temp);%ChowLiu(RHO_temp);
        active_adjacency(index,index)=subgraph;
    end
end

if(n_unconfid_tri==0)
    RHO=X.'*X/size(X,1);
    RHO(logical(eye(p)))=-inf(1,p);
    active_adjacency=ChowLiu(RHO);
end
