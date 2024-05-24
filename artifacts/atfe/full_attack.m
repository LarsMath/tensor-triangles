function Antisymmetrize(F, n, phi)
    return [Matrix(F, n, n, [[phi[i][j][k] + phi[j][k][i] + phi[k][i][j] - phi[i][k][j] - phi[k][j][i] - phi[j][i][k] : i in [1..#phi]] : j in [1..#phi]]) : k in [1..#phi]];
end function;

n := 13;
q := 29;

printf "q: %o, n: %o\n", q, n;
F := GF(q);

// =================================== Find an instance with triangle =======================================

R1<[vars]> := PolynomialRing(F, 3*(n-3), "grevlex");

u := [1,0,0] cat vars[1..n-3];
v := [0,1,0] cat vars[n-2..2*n-6];
w := [0,0,1] cat vars[2*n-5..3*n-9];


// We have to make sure that the orbit of phi has a triangle and that for both phi and psi the triangle does not intersect with the x1=x2=x3=0 hyperplane
// Note that if we really want to break an instance we could lose the not intersecting with x1=x2=x3=0 hyperplane constraint
count := 0;
while true do
    count +:=1;
    phi := Antisymmetrize(F, n, [[[Random(F) : _ in [1..n]] : _ in [1..n]] : _ in [1..n]]);

    system := [&+[phi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]];
    system cat:=[&+[phi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]];
    system cat:=[&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]];

    if VarietySize(Ideal(system)) eq 1 then
        A := Random(GL(n, F));
        A := Matrix(F, n, n, [[A[i][j] : i in [1..n]] : j in [1..n]]);
        psi := [&+[A[i][j] * (Transpose(A) * phi[i] * A) : i in [1..n]] : j in [1..n]];
        system := [&+[psi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]];
        system cat:=[&+[psi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]];
        system cat:=[&+[psi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]];
        if VarietySize(Ideal(system)) eq 1 then
            break;
        end if;
    end if;
end while;

sol := Eltseq(A) cat Eltseq(A^(-1));
print "Found an instance with a triangle in", count, "tries";

start_time := Cputime();

// ================================= Find both triangles ==================================================

system := [&+[phi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]];
system cat:=[&+[phi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]];
system cat:=[&+[phi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]];

phi_triangle := Variety(Ideal(system))[1];

u_phi := [1, 0, 0] cat [phi_triangle[i] : i in [1..n-3]];
v_phi := [0, 1, 0] cat [phi_triangle[i] : i in [n-2..2*n-6]];
w_phi := [0, 0, 1] cat [phi_triangle[i] : i in [2*n-5..3*n-9]];

system := [&+[psi[i][j][l] * v[i] * w[j] : i in [1..n], j in [1..n]] : l in [1..n]];
system cat:=[&+[psi[i][j][l] * w[i] * u[j] : i in [1..n], j in [1..n]] : l in [1..n]];
system cat:=[&+[psi[i][j][l] * u[i] * v[j] : i in [1..n], j in [1..n]] : l in [1..n]];

psi_triangle := Variety(Ideal(system))[1];

u_psi := [1, 0, 0] cat [psi_triangle[i] : i in [1..n-3]];
v_psi := [0, 1, 0] cat [psi_triangle[i] : i in [n-2..2*n-6]];
w_psi := [0, 0, 1] cat [psi_triangle[i] : i in [2*n-5..3*n-9]];

// ===============================Transform instances===================================

while true do
    temp := Matrix(F, n, n, [u_phi, v_phi, w_phi] cat [[Random(F) : _ in [1..n]] : _ in [1..n-3]]);
    if Rank(temp) eq n then
        A_phi := Transpose(temp);
        break;
    end if;
end while;

phi_prime := [&+[A_phi[i][j] * (Transpose(A_phi) * phi[i] * A_phi) : i in [1..n]] : j in [1..n]];

while true do
    temp := Matrix(F, n, n, [u_psi, v_psi, w_psi] cat [[Random(F) : _ in [1..n]] : _ in [1..n-3]]);
    if Rank(temp) eq n then
        A_psi := Transpose(temp);
        break;
    end if;
end while;

psi_prime := [&+[A_psi[i][j] * (Transpose(A_psi) * psi[i] * A_psi) : i in [1..n]] : j in [1..n]];

// ========================= Post collision =================================

wtA12 := 3;
wtA22 := 3;

R2<[x]> := PolynomialRing(F, [1 : _ in [1..9]] cat [wtA12 : _ in [1..(n-3)*3]] cat [wtA22 : _ in [1..(n-3)^2]] cat [1 : _ in [1..9]] cat [wtA12 : _ in [1..(n-3)*3]] cat [wtA22 : _ in [1..(n-3)^2]]);

// Note: Bx representes the inverse of A
Ax := Matrix([(i le 3 select [x[j+3*(i-1)] : j in [1..3]] else [0,0,0]) cat [x[j + (n-3)*(i-1) + 9] : j in [1..n-3]] : i in [1..n]]);
Bx := Matrix([(i le 3 select [x[j+3*(i-1) + 9 + n*(n-3)] : j in [1..3]] else [0,0,0]) cat [x[j + (n-3)*(i-1) + 9 + n*(n-3) + 9] : j in [1..n-3]] : i in [1..n]]);

// ================================================================

In := Matrix(Identity(GL(n,F)));

system_AAB := &cat[Eltseq(Transpose(Ax) * phi_prime[j] * Ax - &+[Bx[i][j] * psi_prime[i] : i in [1..n]]) : j in [1..n]];
system_AAB := Reduce(system_AAB);

system_ABB := &cat[Eltseq(Transpose(Ax) * phi_prime[j] - &+[Bx[i][j] * (psi_prime[i] * Bx) : i in [1..n]]) : j in [1..n]];
system_ABB := Reduce(system_ABB);

system_AB := Eltseq(Ax * Bx - In) cat Eltseq(Bx * Ax - In);
system_AB := Reduce(system_AB);

system := system_AAB cat system_ABB cat system_AB;

s := Variety(Ideal(GroebnerBasis(system, 6)))[1];

As := Matrix([(i le 3 select [s[j+3*(i-1)] : j in [1..3]] else [0,0,0]) cat [s[j + (n-3)*(i-1) + 9] : j in [1..n-3]] : i in [1..n]]);
Bs := Matrix([(i le 3 select [s[j+3*(i-1) + 9 + n*(n-3)] : j in [1..3]] else [0,0,0]) cat [s[j + (n-3)*(i-1) + 9 + n*(n-3) + 9] : j in [1..n-3]] : i in [1..n]]);

print "Found:", A_phi * As * A_psi^(-1);
print "Actual solution:", A;
print "Time since instance found:", Cputime() - start_time;

exit;