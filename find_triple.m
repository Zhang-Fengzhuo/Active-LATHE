function triple=find_triple(A)
%find_triple finds all the triples (i,j,k) that (i,j) and (j,k) are in the
%edge set
%Parameters:
%A: the adjacency matrix of the graph
%
%Output:
%triple: the triples (i,j,k) that (i,j) and (j,k) are in the edge set of the graph
%
%Fengzhuo Zhang, Oct 2021, NUS
triple=[];
[p,~]=size(A);
for i=1:p
    neighbors=find(A(i,:));
    len=length(neighbors);
    if (len==1)
        continue;
    else
        for j=1:len-1
            for k=j+1:len
                triple=[triple;neighbors(j) i neighbors(k)];
            end
        end
    end
end
