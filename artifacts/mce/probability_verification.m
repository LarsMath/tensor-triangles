experiments := 1000;

function CountTriangles(n, m, k, q, experiments)
    F := GF(q);

    R<[vars]> := PolynomialRing(F, (n-1) + (m-1) + (k-1), "grevlex");
    u := [1] cat vars[1..(n-1)];
    v := [1] cat vars[(n-1)+1..(n-1)+(m-1)];
    w := [1] cat vars[(n-1)+(m-1)+1..(n-1)+(m-1)+(k-1)];

    unique := 0;
    non_unique := 0;
    for _ in [1..experiments] do

        // Generate
        phi := [[[Random(F) : _ in [1..k]] : _ in [1..m]] : _ in [1..n]];

        // Build system
        system :=   [&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..m]] : l in [1..k]];
        system cat:=[&+[phi[i][j][l] * v[j] * w[l] : l in [1..k], j in [1..m]] : i in [1..n]];
        system cat:=[&+[phi[i][j][l] * w[l] * u[i] : i in [1..n], l in [1..k]] : j in [1..m]];

        solutions := VarietySize(Ideal(system));
        if solutions eq 1 then
            unique +:= 1;
        elif solutions gt 1 then
            non_unique +:= 1;
        end if;
    end for;
    return unique, non_unique;
end function;

for n in [4..5] do
    for q in [13,31,251] do
        m := n;
        k := n;
        printf "q: %o, n: %o, m: %o, k: %o, experiments: %o\n", q, n, m, k, experiments;
        unique, non_unique := CountTriangles(n, m, k, q, experiments);
        printf "total:%o, unique: %o, non-unique: %o\n", unique + non_unique, unique, non_unique;
    end for;
end for;



