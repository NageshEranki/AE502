function out=h(r20,v20,t,R)

    global data

    rs=zeros(size(data,1)-1,3);
    
    for i=2:size(data,1)
        r=curtis3_4(r20,v20,t(i)-t(2));
        rs(i-1,:)=r;
    end
    out = r2radec(rs-R(2:end,:)); out=out(:);
end

