function TrilinearXLComplexity(hilbert_series, monomial_series, precision)
	least_monomials := 0;
	dsolv := [];
	for dx in [0..precision-1] do
		for dy in [0..precision-1] do
			for dz in [0..precision-1] do
				if Coefficient(Coefficient(Coefficient(hilbert_series, dz), dy), dx) le 0 then
					m := Coefficient(Coefficient(Coefficient(monomial_series, dz), dy), dx);
					if least_monomials eq 0 or m lt least_monomials then
						least_monomials := m;
						dsolv := [<dx, dy, dz>];
                    elif m eq least_monomials then
                        dsolv cat:=[<dx, dy, dz>];
                    end if;
				end if;
			end for;
		end for;
	end for;
	return least_monomials, dsolv;
end function;