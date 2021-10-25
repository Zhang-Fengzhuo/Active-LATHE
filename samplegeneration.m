function X=samplegeneration(adjmat,theta,n)
%samplegeneration generates samples of a homogeneous tree-structured Ising model
%Parameters:
%adjmat: adjacency matrix of the underlying tree graph
%theta: the parameter of the homogeneous Ising tree model
%n: the number of samples to generate
%
%Output:
%X: n*p data matrix
%
%Fengzhuo Zhang, Oct 2021, NUS
root=1;
[p,~]=size(adjmat);
p_cross=theta;%exp(-theta)/(exp(theta)+exp(-theta));
X=zeros(n,p);
for i=1:n
    x_r=(rand>0.5)*2-1;
    X(i,root)=x_r;
    child=find(adjmat(root,:)>0);
    pa=root*ones(1,length(child));
    while(~isempty(child))
        pa_next=[];
        child_next=[];
        for j=1:length(child)
            x=X(i,pa(j))*((rand>p_cross)*2-1);
            X(i,child(j))=x;
            temp_c=find(adjmat(child(j),:)>0);
            temp_c=temp_c(find(temp_c~=pa(j)));
            child_next=[child_next temp_c];
            pa_next=[pa_next child(j)*ones(1,length(temp_c))];
        end
        child=child_next;
        pa=pa_next;
    end
end