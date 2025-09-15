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
using CodecXz
using ProgressBars

@assert length(ARGS) == 1 "Usage: parse-triangulations.batch.jl folder"
folder, = ARGS
files = readdir(folder; join=true) |> filter(contains(r"\d.topcom"))

for dat_file in ProgressBar(files)
  name    = rsplit(dat_file, "."; limit=2)[1]
  xz_file = name * ".out.xz"
  outfile = name * ".mrdi"

  entire = read(dat_file, String)
  m = match(r"(\[\D*\[.*\]\D*\])\D*(\[\D*\[.*\]\D*\])"s, entire)
  if !isnothing(m)
    points = m.captures[1]
  else
    points = entire
  end
  points = matrix(QQ, Meta.eval(Meta.parse(points)))
  try
    global group = m.captures[2]
    global group = Polymake.to_one_based_indexing(convert(Vector{Vector{Int}}, Meta.parse(group)|>Meta.eval))
    global group = permutation_group(nrows(points), [perm(x) for x in group])
  catch e
    global group = permutation_group(nrows(points), PermGroupElem[])
  end

  triangs = IncidenceMatrix[]
  global i = 1
  io = XzDecompressorStream(open(xz_file,"r");memlimit=81920000)
  while !eof(io)
      line = readline(io)
      m = match(r"{{.*}}", line)
      triang = replace(m.match, "{" => "[")
      triang = replace(triang, "}" => "]")
      triang = convert(Vector{Vector{Int}}, JSON.parse(triang))
      triang = Polymake.to_one_based_indexing(triang)
      push!(triangs, IncidenceMatrix(triang))
      global i += 1
  end
  close(io)

  println("Found $(i-1) triangulations in $xz_file")
  tuple = (points=points, group=group, triangulations=triangs)
  out = XzCompressorStream(open(outfile, "w+"); memlimit=81920000)
    save(out, tuple)
  close(out)
end

