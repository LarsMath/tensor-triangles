load "../tools.m";
load "./find_triangles_algorithm.m";

function TestTriangles(n, q, experiments)
    F := GF(q);

    unique := 0;
    non_unique := 0;
    for _ in [1..experiments] do

        // Generate random ATF
        phi := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);

        solutions := #FindTriangles(phi);
        if solutions eq 1 then
            unique +:= 1;
        elif solutions gt 1 then
            non_unique +:= 1;
        end if;
    end for;
    return unique, non_unique;
end function;
    
experiments := 1000;
n := 9;
q := 29;

printf "q: %o, n: %o, experiments: %o\n", q, n, experiments;
unique, non_unique := TestTriangles(n, q, experiments);
printf "experiments: %o, unique: %o, non-unique: %o\n", experiments, unique, non_unique;

exit;