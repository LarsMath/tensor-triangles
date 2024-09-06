load "./find_triangles_algorithm.m";

function TestTriangles(n, m, k, q, experiments)
    F := GF(q);

    unique := 0;
    non_unique := 0;
    for _ in [1..experiments] do

        // Generate random tensor
        C := [[[Random(F) : _ in [1..k]] : _ in [1..m]] : _ in [1..n]];

        solutions := #FindTriangles(C);
        if solutions eq 1 then
            unique +:= 1;
        elif solutions gt 1 then
            non_unique +:= 1;
        end if;
    end for;
    return unique, non_unique;
end function;

experiments := 1000;

q := 31;
n := 5;
m := 5;
k := 5;

printf "q: %o, n: %o, m: %o, k: %o, experiments: %o\n", q, n, m, k, experiments;
unique, non_unique := TestTriangles(n, m, k, q, experiments);
printf "experiments: %o, unique: %o, non-unique: %o\n", experiments, unique, non_unique;

exit;
