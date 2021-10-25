function CL_adjacency=ChowLiu(mutual_info)
%ChowLiu estimates the structure of model by Chow-Liu algorithm
%Parameters:
%mutual_info: mutual information matrix (It can be the matrix of correlation in Ising model with zero-field, see Tandon, Tan, and Zhu 2020 and Bresler and Karzand 2020)
%
%Output:
%CL_adjacency: the adjacency matrix estimate
%
%Fengzhuo Zhang, Oct 2021, NUS
[p,~]=size(mutual_info);
CL_adjacency=zeros(p,p);
V=[1:p];
U=[1];
for i=1:p-1
    W=setdiff(V,U);
    sub=mutual_info(U,W);
    [m,n]=find(sub==max(max(sub)));
    U=[U W(n(1))];
    CL_adjacency(U(m(1)),W(n(1)))=1;
end
CL_adjacency=CL_adjacency+CL_adjacency.';
    
