function Dy=diffy(D)
[m,n]=size(D);
Dy=zeros(m,m);
        for j=1:m-1
            Dy(j,j)=1;
            Dy(j,j+1)=-1;
        end
        Dy(m,m)=1;
        Dy(m,1)=-1;
end
