load "../tools.m";
load "./post_collision_algorithm.m";

q := 2^32 - 5;
n := 13;

printf "Post-collision (n: %o, q: %o)\n", n, q;

F := GF(q);

// Generate problem
while true do
    A := Matrix(F, n, n, [[i gt 3 and j le 3 select 0 else Random(F) : j in [1..n]] : i in [1..n]]);
    if Rank(A) eq n then break; end if;
end while;

phi := Antisymmetrize([[[#[a : a in [i,j,k] | a le 3] ge 2 select 0 else Random(F) : k in [1..n]] : j in [1..n]] : i in [1..n]]);
psi := ApplyGroupAction3TI(A, A, A, phi, F);

B := A^(-1);
solution := [A[i][j] : j in [1..3], i in [1..3]] cat [A[i][j] : j in [4..n], i in [1..n]] cat [B[i][j] : j in [1..3], i in [1..3]] cat [B[i][j] : j in [4..n], i in [1..n]];

found := PostCollision(phi, psi : verbosity:=1);

if #found eq 0 then
    print "This should not happen";
else
    print "Found constructed solution:", solution in [[x : x in f] : f in found];
end if;

exit;