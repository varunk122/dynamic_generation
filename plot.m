T = readtable("values.csv");
[x y z] = sphere;

figure;
hold on;
for i= 1:height(T)
    x_ = T.Var1(i);
    y_ = T.Var2(i);
    z_ = T.Var3(i);
    r = T.Var4(i);
    
    surf(x*r+x_,y*r+y_,z*r+z_);
end

grid on;
daspect([1 1 1])
view(30,10)