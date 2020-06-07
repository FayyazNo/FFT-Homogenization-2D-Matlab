function multpCx =gg(Ci, Cm, RVE, ep)

X=size(ep,1);
Y=size(ep,2);
multpCx=zeros(size(ep));
for i=1:X
    for j=1:Y
      if RVE(i,j)==0 
         
         multpCx(i,j,:)=Cm*reshape(ep(i,j,:),3,1);

      else
        
         multpCx(i,j,:)=Ci*reshape(ep(i,j,:),3,1);
      end
    
    end
    
end
