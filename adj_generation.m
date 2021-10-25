function [edge,Adj]=adj_generation(name,p)
%adj_generation generates the adjacency matrix of the desired model
%Parameters:
%name: the name of the desired model
%p: the number of nodes
%
%Output:
%Adj: the adjacency matrix of the desired model
%edge: p*2 matrix of the edges in the graph
%
%Fengzhuo Zhang, Oct 2021, NUS

edge=[];
Adj=zeros(p,p);

switch name
    case 'chain'
        for i=1:p-1
            edge=[edge;i i+1;];
        end
        [N_e,~]=size(edge);
        Adj=zeros(p,p);
        for i=1:N_e
            Adj(edge(i,1),edge(i,2))=1;
            Adj(edge(i,2),edge(i,1))=1;
        end
        
    case 'hmm'
        if(mod(p,2)==0)
            n_layer1=round(p/2)+1;
            for i=1:n_layer1-1
                edge=[edge; i i+1];
                Adj(i,i+1)=1;
                Adj(i+1,i)=1;
            end
            index1=[2:(n_layer1-1)];
            index2=[(n_layer1+1):p];
            new_edge=[index1.' index2.'];
            edge=[edge;new_edge];
            for i=1:length(index1)
                Adj(index1(i),index2(i))=1;
                Adj(index2(i),index1(i))=1;
            end
                
        end
        
    case 'binarytree'
        L=round(log2(p+1));
        for i=2:L
            pa=[2^(i-2):(2^(i-1)-1)];
            index_node=2^(i-1):min([p 2^i-1]);
            num_node=length(index_node);
            for j=1:num_node
                edge=[edge; index_node(j) pa(ceil(j/2))];
                Adj(index_node(j),pa(ceil(j/2)))=1;
                Adj(pa(ceil(j/2)),index_node(j))=1;
            end
        end
    otherwise
        disp('No such model.')
end
        