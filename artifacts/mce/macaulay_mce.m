load "../tools.m";

q := 29;

n := 7;
m := 7;
k := 7;

dx := 2;
dy := 2;
dz := 1;

F := GF(q);

printf "n, m, k: %o, %o, %o\tq: %o\tdx, dy, dz: %o, %o, %o\n", n, m, k, q, dx, dy, dz;


R<[vars]> := PolynomialRing(F, (n-1) + (m-1) + (k-1), "grevlex");
xs := vars[1..(n-1)];
ys := vars[(n-1)+1..(n-1)+(m-1)];
zs := vars[(n-1)+(m-1)+1..(n-1)+(m-1)+(k-1)];

phi := [[[Random(F) : i in [1..k]] : j in [1..m]] : l in [1..n]];

system_xy := Reduce([&+[phi[i][j][l] * xs[i] * ys[j] : i in [1..n-1], j in [1..m-1]] : l in [1..k]]);
system_yz := Reduce([&+[phi[i][j][l] * ys[j] * zs[l] : l in [1..k-1], j in [1..m-1]] : i in [1..n]]);
system_zx := Reduce([&+[phi[i][j][l] * zs[l] * xs[i] : i in [1..n-1], l in [1..k-1]] : j in [1..m]]);

mxs := GenerateMonomials(xs, dx);
mys := GenerateMonomials(ys, dy);
mzs := GenerateMonomials(zs, dz);

system :=   [f * x * y * z : f in system_xy, x in mxs[dx+1], y in mys[dy+1], z in mzs[dz+2]];
system cat:=[f * x * y * z : f in system_yz, x in mxs[dx+2], y in mys[dy+1], z in mzs[dz+1]];
system cat:=[f * x * y * z : f in system_zx, x in mxs[dx+1], y in mys[dy+2], z in mzs[dz+1]];

ms := Reverse(Sort([x*y*z : x in mxs[dx+2], y in mys[dy+2], z in mzs[dz+2]]));
system := Reverse(Sort(system));

printf "Building and reducing %o x %o matrix\n", #system, #ms;
M := Matrix(F, #system, #ms, [[MonomialCoefficient(f, m) : m in ms] : f in system]);

print "Syzygies:", #system - Rank(M);

exit;
nb