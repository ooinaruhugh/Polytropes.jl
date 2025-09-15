###############################################################################
###############################################################################
# parsing-triangulations.jl
#
# This file parses output of TOPCOM or mptopcom into Julia and OSCAR.
# Additionally it also reads the input file from TOPCOM or mptopcom to make
# sure that one works with the same points and group as TOPCOM or mptopcom.
# Furthermore, this script converts to the appropriate OSCAR datatypes, e.g.
# the group becomes an OSCAR group and the triangulations have the type
# `SubdivisionOfPoints`.
#
# Usage:
# julia parsing-triangulations.jl DAT_FILE XZ_FILE OUT_FILE
#
# Since output from TOPCOM or mptopcom can be very large it is assumed to be
# given in a compressed form.
#
# Examples:
# julia parsing-triangulations.jl D4xD2.dat points2triangs.out.xz points2triangs.mrdi
# julia parsing-triangulations.jl D4xD2.dat mptopcom1.out.xz mptopcom1.mrdi
# julia parsing-triangulations.jl D4xD2.dat mptopcom.out.xz mptopcom.mrdi
#
###############################################################################
###############################################################################
using JSON
using Oscar

@assert length(ARGS) == 3 "Usage: parse-triangulations.jl DAT_FILE XZ_FILE OUT_FILE"

include("parse-triangulations.jl")
parse_triangulation_to_mrdi(ARGS...)
