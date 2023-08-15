variable runnerCutoff    equal  6.352  # largest symmetry function cutoff (Angstrom)
variable runnerDir    string "../nnp/"

mass 3 47.86
mass 2 16.00
mass 1 1.00

pair_style nnp dir ${runnerDir} showew no showewsum 10 resetew yes maxew 500 cflength 1.8897261328 cfenergy 0.0015936014
pair_coeff * * ${runnerCutoff}        # set up pair style coefficients
