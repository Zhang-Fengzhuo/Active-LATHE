function [CP,n_CP]=find_commedge_component(triple)
%find_commedge_component finds the connected components of the nodes and
%edges in variable 'triple'
%Parameters:
%triple: a set of triples (i,j,k), where i,j,k are nodes and i-j and j-k
%are edges
%
%Output:
%CP: the connected components on the nodes and edges in variable 'triple'
%n_CP: the number of the connected components
%
%Fengzhuo Zhang, Oct 2021, NUS
[n_triple,~]=size(triple);
tempCP=cell(n_triple,1);
n_CP=0;
for i=1:n_triple
    flag=0;
    for j=1:n_CP
        if(length(intersect(triple(i,:), tempCP{j,1}))>=2)
            tempCP{j,1}=union(tempCP{j,1},triple(i,:));
            flag=1;
            break;
        end
    end
    if(flag==0)
        tempCP{n_CP+1,1}=triple(i,:);
        n_CP=n_CP+1;
    end
end
CP=cell(n_CP,1);
for i=1:n_CP
    CP{i,1}=tempCP{i,1};
end
