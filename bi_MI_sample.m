function MI=bi_MI_sample(data)
%bi_MI_sample calculates the mutual information between two binary
%variables from their samples
%Parameters:
%data is a n*2 matrix with entry -1,1
%
%Output:
%MI: the mutual information between nodes
%
%Fengzhuo Zhang, Oct 2021, NUS
[len,~]=size(data);
p1=sum(data(:,1)==1)/len;
p2=sum(data(:,2)==1)/len;
p11=sum((data(:,1)==1).*(data(:,2)==1))/len;
p01=sum((data(:,1)==-1).*(data(:,2)==1))/len;
p10=sum((data(:,1)==1).*(data(:,2)==-1))/len;
p00=sum((data(:,1)==-1).*(data(:,2)==-1))/len;

%MI=p11+p00;

MI=0;
if(p11>0)
    MI=MI+p11*log(p11/(p1*p2));
end
if(p01>0)
    MI=MI+p01*log(p01/((1-p1)*p2));
end
if(p10>0)
    MI=MI+p10*log(p10/(p1*(1-p2)));
end
if(p00>0)
    MI=MI+p00*log(p00/((1-p1)*(1-p2)));
end
%}

