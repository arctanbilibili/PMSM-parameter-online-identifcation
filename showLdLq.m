file=fullfile("./", "xk6.txt_MARS.txt");
x6kdata = load(file);
ss =  size(x6kdata);


x6begin = 1;
x6end   = ss(1);

%多少列
x6col = ss(2);

rng(114);
nois01 = 0;
nois02 = 0;

figure(2)

for i=1:x6col
subplot(x6col,1,i)
hold on
plot(x6kdata(:,i));
grid on;
hold off
end




