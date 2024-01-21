using Oscar
using Polytropes

G = complete_dag(3)
I = toric_ideal(G)
GF,gbs = groebner_fan(I; return_groebner_bases=true)
# GF,gbs,ords,initial_ideals = groebner_fan(I; return_groebner_bases=true,return_orderings=true,return_initial_ideals=true)