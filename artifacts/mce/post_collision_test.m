load "../tools.m";
load "./post_collision_algorithm.m";

q := 4093;
n := 14;
m := 14;
k := 14;

printf "Post-collision (n: %o, m: %o, k: %o, q: %o)\n", n, m, k, q;

F := GF(q);

// Construct problem with triangles in standard position 

C := [[[Random(F) : _ in [1..k]] : _ in [1..m]] : _ in [1..n]];

for i in [1..n] do C[i][1][1] := 0; end for;
for j in [1..m] do C[1][j][1] := 0; end for;
for l in [1..k] do C[1][1][l] := 0; end for;

e1A := Matrix(F, n, 1, [1] cat [0 : _ in [1..n-1]]);
e1B := Matrix(F, n, 1, [1] cat [0 : _ in [1..m-1]]);
e1T := Matrix(F, n, 1, [1] cat [0 : _ in [1..k-1]]);

while true do
    A := HorizontalJoin(Random(F) * e1A, RandomMatrix(F, n, n-1));
    B := HorizontalJoin(e1B, RandomMatrix(F, m, m-1));
    T := HorizontalJoin(e1T, RandomMatrix(F, k, k-1));
    if Rank(A) eq n and Rank(B) eq m and Rank(T) eq k then break; end if;
end while;

D := ApplyGroupAction3TI(A, B, T, C, F);

solution := Eltseq(ColumnSubmatrix(A, 2, n-1)) cat Eltseq(ColumnSubmatrix(A^(-1), 2, n-1)) cat Eltseq(ColumnSubmatrix(B, 2, m-1)) cat Eltseq(ColumnSubmatrix(B^(-1), 2, m-1)) cat Eltseq(ColumnSubmatrix(T, 2, k-1)) cat Eltseq(ColumnSubmatrix(T^(-1), 2, k-1)) cat [A[1][1], 1/A[1][1]];

found := PostCollision(C, D : verbosity:=1);

if #found eq 0 then
    print "This should not happen";
else
    print "Found constructed solution:", solution in [[x : x in f] : f in found];
end if;

exit;
