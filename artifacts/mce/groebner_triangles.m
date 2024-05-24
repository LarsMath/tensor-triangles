q := 13;
n := 9;
m := 8;
k := 8;

F := GF(q);

R<[vars]> := PolynomialRing(F, (n-1) + (m-1) + (k-1), "grevlex");
u := [1] cat vars[1..(n-1)];
v := [1] cat vars[(n-1)+1..(n-1)+(m-1)];
w := [1] cat vars[(n-1)+(m-1)+1..(n-1)+(m-1)+(k-1)];

phi := [[[Random(F) : i in [1..k]] : j in [1..m]] : l in [1..n]];

system :=   [&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..m]] : l in [1..k]];
system cat:=[&+[phi[i][j][l] * v[j] * w[l] : l in [1..k], j in [1..m]] : i in [1..n]];
system cat:=[&+[phi[i][j][l] * w[l] * u[i] : i in [1..n], l in [1..k]] : j in [1..m]];

SetVerbose("Groebner", 1);
GB := GroebnerBasis(system);

exit;


