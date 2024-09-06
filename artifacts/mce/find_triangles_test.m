load "./find_triangles_algorithm.m";

q := 4093;
n := 7;
m := 7;
k := 7;

printf "Find triangles (n: %o, m: %o, k: %o, q: %o)\n", n, m, k, q;

F := GF(q);

C := [[[Random(F) : _ in [1..k]] : _ in [1..m]] : _ in [1..n]];

FindTriangles(C : verbosity := 1);

exit;
