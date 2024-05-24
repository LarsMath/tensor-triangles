load "../tools.m";

q := 31;
n := 8;
m := n;
k := n;
print "n:", n, "m:", m, "k:", k;

F := GF(q);
R<[x]> := PolynomialRing(F, 2 * (1 + n*(n-1) + m*(m-1) + k*(k-1)), "grevlex");

// Construct problem with triangles in standard position 

phi := [[[Random(F) : i in [1..n]] : j in [1..n]] : k in [1..n]];
for i in [1..n] do
    phi[i][1][1] := 0;
    phi[1][i][1] := 0;
    phi[1][1][i] := 0;
end for;

A_sol := HorizontalJoin(Matrix(R, n, 1, [[i eq 1 select Random(F) else 0] : i in [1..n]]), Matrix(R, n, n-1, [[Random(F) : _ in [1..n-1]] : _ in [1..n]]));
B_sol := HorizontalJoin(Matrix(R, m, 1, [[i eq 1 select 1 else 0] : i in [1..m]]), Matrix(R, m, m-1, [[Random(F) : _ in [1..m-1]] : _ in [1..m]]));
T_sol := HorizontalJoin(Matrix(R, k, 1, [[i eq 1 select 1 else 0] : i in [1..k]]), Matrix(R, k, k-1, [[Random(F) : _ in [1..k-1]] : _ in [1..k]]));
while Rank(A_sol) ne n or Rank(B_sol) ne m or Rank(T_sol) ne k or A_sol[1][1] eq 0 do
    A_sol := HorizontalJoin(Matrix(R, n, 1, [[i eq 1 select Random(F) else 0] : i in [1..n]]), Matrix(R, n, n-1, [[Random(F) : _ in [1..n-1]] : _ in [1..n]]));
    B_sol := HorizontalJoin(Matrix(R, m, 1, [[i eq 1 select 1 else 0] : i in [1..m]]), Matrix(R, m, m-1, [[Random(F) : _ in [1..m-1]] : _ in [1..m]]));
    T_sol := HorizontalJoin(Matrix(R, k, 1, [[i eq 1 select 1 else 0] : i in [1..k]]), Matrix(R, k, k-1, [[Random(F) : _ in [1..k-1]] : _ in [1..k]]));
end while;

psi := ApplyGroupAction3TI(A_sol, B_sol, T_sol, phi, R);

sol := Eltseq(ColumnSubmatrix(A_sol, 2, n-1)) cat Eltseq(ColumnSubmatrix(A_sol^(-1), 2, n-1)) cat Eltseq(ColumnSubmatrix(B_sol, 2, m-1)) cat Eltseq(ColumnSubmatrix(B_sol^(-1), 2, m-1)) cat Eltseq(ColumnSubmatrix(T_sol, 2, k-1)) cat Eltseq(ColumnSubmatrix(T_sol^(-1), 2, k-1)) cat [A_sol[1][1], 1/A_sol[1][1]];

// ========================================================

start_n_time := Cputime();

// Build variable matrices
A := HorizontalJoin(Matrix(R, n, 1, [[i eq 1 select x[2*(n*(n-1) + m*(m-1) + k*(k-1))+1] else 0] : i in [1..n]]), Matrix(R, n, n-1, x[1..n*(n-1)]));
Ai := HorizontalJoin(Matrix(R, n, 1, [[i eq 1 select x[2*(n*(n-1) + m*(m-1) + k*(k-1))+2] else 0] : i in [1..n]]), Matrix(R, n, n-1, x[n*(n-1)+1..2*n*(n-1)]));
B := HorizontalJoin(Matrix(R, m, 1, [[i eq 1 select 1 else 0] : i in [1..m]]), Matrix(R, m, m-1, x[2*n*(n-1)+1..2*n*(n-1)+m*(m-1)]));
Bi := HorizontalJoin(Matrix(R, m, 1, [[i eq 1 select 1 else 0] : i in [1..m]]), Matrix(R, m, m-1, x[2*n*(n-1)+m*(m-1)+1..2*n*(n-1)+2*m*(m-1)]));
T := HorizontalJoin(Matrix(R, k, 1, [[i eq 1 select 1 else 0] : i in [1..k]]), Matrix(R, k, k-1, x[2*n*(n-1)+2*m*(m-1)+1..2*n*(n-1)+2*m*(m-1)+k*(k-1)]));
Ti := HorizontalJoin(Matrix(R, k, 1, [[i eq 1 select 1 else 0] : i in [1..k]]), Matrix(R, k, k-1, x[2*n*(n-1)+2*m*(m-1)+k*(k-1)+1..2*n*(n-1)+2*m*(m-1)+2*k*(k-1)]));


//=========================== Build system ========================================

In := Matrix(R, n, n, [[i eq j select 1 else 0 : i in [1..n]] : j in [1..n]]);
Im := Matrix(R, m, m, [[i eq j select 1 else 0 : i in [1..m]] : j in [1..m]]);
Ik := Matrix(R, k, k, [[i eq j select 1 else 0 : i in [1..k]] : j in [1..k]]);

phi_AB := ApplyGroupAction3TI(A, B, Ik, phi, R);
phi_AT := ApplyGroupAction3TI(A, Im, T, phi, R);
phi_BT := ApplyGroupAction3TI(In, B, T, phi, R);
phi_A := ApplyGroupAction3TI(A, Im, Ik, phi, R);
phi_B := ApplyGroupAction3TI(In,  B, Ik, phi, R);
phi_T := ApplyGroupAction3TI(In, Im, T, phi, R);

psi_AB := ApplyGroupAction3TI(Ai, Bi, Ik, psi, R);
psi_AT := ApplyGroupAction3TI(Ai, Im, Ti, psi, R);
psi_BT := ApplyGroupAction3TI(In, Bi, Ti, psi, R);
psi_A := ApplyGroupAction3TI(Ai, Im, Ik, psi, R);
psi_B := ApplyGroupAction3TI(In, Bi, Ik, psi, R);
psi_T := ApplyGroupAction3TI(In, Im, Ti, psi, R);

AAi := Matrix(A) * Matrix(Ai);
BBi := Matrix(B) * Matrix(Bi);
TTi := Matrix(T) * Matrix(Ti);
AiA := Matrix(Ai) * Matrix(A);
BiB := Matrix(Bi) * Matrix(B);
TiT := Matrix(Ti) * Matrix(T);

system := [];
system := system cat [phi_A[i][j][ii] - psi_BT[i][j][ii] : i in [1..n], j in [1..m], ii in [1..k]];
system := system cat [phi_B[i][j][ii] - psi_AT[i][j][ii] : i in [1..n], j in [1..m], ii in [1..k]];
system := system cat [phi_T[i][j][ii] - psi_AB[i][j][ii] : i in [1..n], j in [1..m], ii in [1..k]];
system := system cat [phi_AB[i][j][ii] - psi_T[i][j][ii] : i in [1..n], j in [1..m], ii in [1..k]];
system := system cat [phi_AT[i][j][ii] - psi_B[i][j][ii] : i in [1..n], j in [1..m], ii in [1..k]];
system := system cat [phi_BT[i][j][ii] - psi_A[i][j][ii] : i in [1..n], j in [1..m], ii in [1..k]];

system := system cat [AAi[i][j] - In[i][j] : i in [1..n], j in [1..n]];
system := system cat [BBi[i][j] - Im[i][j] : i in [1..m], j in [1..m]];
system := system cat [TTi[i][j] - Ik[i][j] : i in [1..k], j in [1..k]];
system := system cat [AiA[i][j] - In[i][j] : i in [1..n], j in [1..n]];
system := system cat [BiB[i][j] - Im[i][j] : i in [1..m], j in [1..m]];
system := system cat [TiT[i][j] - Ik[i][j] : i in [1..k], j in [1..k]];

print "linears", #Reduce([f : f in system | Degree(f) eq 1]);
print "quadratics", #Reduce([f : f in system | Degree(f) eq 2]); // This takes long and is not necessary for computation

start_gb_time := Cputime();
SetVerbose("Groebner", 1);
GB := GroebnerBasis(system);

print "solved, total time:", Cputime() - start_n_time, "GB time:", Cputime() - start_gb_time, "\n";

print "Found:", Variety(Ideal(system));

exit;
