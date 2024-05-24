function Antisymmetrize(F, n, phi)
    return [Matrix(F, n, n, [[phi[i][j][k] + phi[j][k][i] + phi[k][i][j] - phi[i][k][j] - phi[k][j][i] - phi[j][i][k] : i in [1..#phi]] : j in [1..#phi]]) : k in [1..#phi]];
end function;

F := GF(29);
n := 13;

print "n = ", n;

// Generate problem
while true do
    A := Matrix(F, n, n, [[i gt 3 and j le 3 select 0 else Random(F) : j in [1..n]] : i in [1..n]]);
    if Rank(A) eq n then
        B := A^(-1);
        break;
    end if;
end while;

phi := Antisymmetrize(F, n, [[[#[a : a in [i,j,k] | a le 3] ge 2 select 0 else Random(F) : k in [1..n]] : j in [1..n]] : i in [1..n]]);
psi := [&+[A[i][j] * (Transpose(A) * phi[i] * A) : i in [1..n]] : j in [1..n]];

In := Matrix(Identity(GL(n,F)));

// ========================================================================

wtA12 := 3;
wtA22 := 3;

R<[x]> := PolynomialRing(F, [1 : _ in [1..9]] cat [wtA12 : _ in [1..(n-3)*3]] cat [wtA22 : _ in [1..(n-3)^2]] cat [1 : _ in [1..9]] cat [wtA12 : _ in [1..(n-3)*3]] cat [wtA22 : _ in [1..(n-3)^2]]);

// Note: Bx representes the inverse of A
Ax := Matrix([(i le 3 select [x[j+3*(i-1)] : j in [1..3]] else [0,0,0]) cat [x[j + (n-3)*(i-1) + 9] : j in [1..n-3]] : i in [1..n]]);
Bx := Matrix([(i le 3 select [x[j+3*(i-1) + 9 + n*(n-3)] : j in [1..3]] else [0,0,0]) cat [x[j + (n-3)*(i-1) + 9 + n*(n-3) + 9] : j in [1..n-3]] : i in [1..n]]);

sol := [A[i][j] : j in [1..3], i in [1..3]] cat [A[i][j] : j in [4..n], i in [1..n]] cat [B[i][j] : j in [1..3], i in [1..3]] cat [B[i][j] : j in [4..n], i in [1..n]];

// ================================================================

system_AAB := &cat[Eltseq(Transpose(Ax) * phi[j] * Ax - &+[Bx[i][j] * psi[i] : i in [1..n]]) : j in [1..n]];
system_AAB := Reduce(system_AAB);
printf "AAB   length: %o\n", {* Degree(f) : f in system_AAB *};

system_ABB := &cat[Eltseq(Transpose(Ax) * phi[j] - &+[Bx[i][j] * (psi[i] * Bx) : i in [1..n]]) : j in [1..n]];
system_ABB := Reduce(system_ABB);
printf "ABB   length: %o\n", {* Degree(f) : f in system_ABB *};

system_AB := Eltseq(Ax * Bx - In) cat Eltseq(Bx * Ax - In);
system_AB := Reduce(system_AB);
printf "AB    length: %o\n", {* Degree(f) : f in system_AB *};

system := system_AAB cat system_ABB cat system_AB;
printf "Total length: %o\n", {* Degree(f) : f in system *};

// ============================================

print "Correct:", [Evaluate(f, sol) : f in system];

// =============================================

SetVerbose("Groebner", 1);
print Variety(Ideal(GroebnerBasis(system, 6)));
exit;