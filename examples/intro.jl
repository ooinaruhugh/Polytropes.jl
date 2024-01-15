using Oscar
using Polytropes

G = complete_dag(3)
I = toric_ideal(G)
GF = groebner_fan(I)