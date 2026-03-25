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

if length(ARGS) != 1 
  println("Usage: parse-triangulations.batch.jl folder")
  exit(1)
end

using JSON
using Oscar
using CodecXz
using ProgressBars

folder, = ARGS
files = readdir(folder; join=true) |> filter(contains(r"\d.topcom"))

include("parse-triangulations.jl")

for dat_file in ProgressBar(files)
  name    = rsplit(dat_file, "."; limit=2)[1]
  xz_file = name * ".out.xz"
  outfile = name * ".mrdi.xz"

  parse_triangulation_to_mrdi(dat_file, xz_file, outfile)
end

