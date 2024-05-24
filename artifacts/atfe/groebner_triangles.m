load "../tools.m";

n := 13;
q := 29;

printf "q: %o, n: %o\n", q, n;

F := GF(q);
R<[vars]> := PolynomialRing(F, 3*(n-3), "grevlex");
u := [1,0,0] cat vars[1..n-3];
v := [0,1,0] cat vars[n-2..2*n-6];
w := [0,0,1] cat vars[2*n-5..3*n-9];

phi := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);

system := [&+[phi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]];
system cat:=[&+[phi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]];
system cat:=[&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]];

SetVerbose("Groebner", 1);
GB := GroebnerBasis(system);

exit;
