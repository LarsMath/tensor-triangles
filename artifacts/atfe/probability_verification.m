load "../tools.m";

experiments := 1000;
n := 9;
q := 251;

function CountTriangles(n,q, experiments)

    F := GF(q);
    R<[x]> := PolynomialRing(F, 3*(n-3), "grevlex");
    v := [1,0,0] cat x[1..n-3];
    w := [0,1,0] cat x[n-2..2*n-6];
    u := [0,0,1] cat x[2*n-5..3*n-9];

    unique := 0;
    non_unique := 0;
    for _ in [1..experiments] do

        // Generate
        phi := Antisymmetrize([[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);

        // Build system
        system :=   [&+[phi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]];
        system cat:=[&+[phi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]];
        system cat:=[&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]];


        solutions := VarietySize(Ideal(system));
        if solutions eq 1 then
            unique +:= 1;
        elif solutions gt 1 then
            non_unique +:= 1;
        end if;
    end for;
    return unique, non_unique;
end function;
    

for n in [7..9] do
    for q in [13,31,251] do
        printf "q: %o, n: %o, experiments: %o\n", q, n, experiments;
        unique, non_unique := CountTriangles(n, q, experiments);
        printf "total:%o, unique: %o, non-unique: %o\n", unique + non_unique, unique, non_unique;
    end for;
end for;

exit;