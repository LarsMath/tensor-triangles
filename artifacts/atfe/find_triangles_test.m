load "../tools.m";
load "./find_triangles_algorithm.m";

n := 13;
q := 2^32 - 5;

printf "Find triangles (n: %o, q: %o)\n", n, q;

F := GF(q);

phi := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);

FindTriangles(phi : verbosity := 1);

exit;
