function unconfid_triple=find_unconfid_triple(triple,RHO,rho)
% find_unconfid_triple finds all unconfident triples.
% Parameters:
% triple: the triples (i,j,k) that (i,j) and (j,k) are in the edge set of the graph
% RHO: correlation matrix where RHO_{i,j} is E[X_{i}X_{j}]
% rho: the estimate of the correlation between adajacent nodes
%
%Output:
% unconfid_triple: triples (i,j,k) that are unconfident
%
%Fengzhuo Zhang, Oct 2021, NUS
unconfid_triple=[];
[len,~]=size(triple);
for i=1:len
    if((RHO(triple(i,1),triple(i,3))>RHO(triple(i,1),triple(i,2))*(11+9*rho)/20)||(RHO(triple(i,1),triple(i,3))>RHO(triple(i,2),triple(i,3))*(11+9*rho)/20))
        unconfid_triple=[unconfid_triple;triple(i,:)];
    end
end
