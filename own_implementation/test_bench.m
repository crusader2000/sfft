n=256;
k=4;
signal=randi([0,n],1,k);
x=zeros(1,n);
x(signal)=100;
plot(1:n,x);