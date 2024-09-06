load "../tools.m";

q := 29;
n := 13;

dx := 2;
dy := 1;
dz := 1;

printf "n: %o\tq: %o\tdx, dy, dz: %o, %o, %o\n", n, q, dx, dy, dz;

F := GF(q);
R<[vars]> := PolynomialRing(F, 3*(n-2), "grevlex");
xs := vars[1..n-2];
ys := vars[n-1..2*n-4];
zs := vars[2*n-3..3*n-6];

v := [xs[1],0,0] cat xs[2..n-2];
w := [0,ys[1],0] cat ys[2..n-2];
u := [0,0,zs[1]] cat zs[2..n-2];

mxs := GenerateMonomials(xs, dx);
mys := GenerateMonomials(ys, dy);
mzs := GenerateMonomials(zs, dz);

ms := Reverse(Sort([x*y*z : x in mxs[dx+2], y in mys[dy+2], z in mzs[dz+2]]));
msi := AssociativeArray(ms);
t := 1; for x in ms do msi[x] := t; t +:= 1; end for;

//===========================

phi := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);

system_xy := Reduce([&+[phi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]]);
system_yz := Reduce([&+[phi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]]);
system_zx := Reduce([&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]]);

system :=   [f * x * y * z : f in system_xy, x in mxs[dx+1], y in mys[dy+1], z in mzs[dz+2]];
system cat:=[f * x * y * z : f in system_yz, x in mxs[dx+2], y in mys[dy+1], z in mzs[dz+1]];
system cat:=[f * x * y * z : f in system_zx, x in mxs[dx+1], y in mys[dy+2], z in mzs[dz+1]];
system := Reverse(Sort(system));

printf "Building and reducing %o x %o matrix\n", #system, #ms;
M := SparseMatrix(F, #system, #ms, [<i, msi[j], MonomialCoefficient(system[i], j)> : j in Monomials(system[i]), i in [1..#system]]);

print "Syzygies:", #system - Rank(M);

exit;
