variable runnerCutoff    equal  6.352  # largest symmetry function cutoff (Angstrom)

mass 3 47.86
mass 2 16.00
mass 1 1.00

#pair_style nnp dir ${nnDir} showew no showewsum 10 resetew yes maxew 100000 cflength 1.8897261328 cfenergy 0.0015936014
#pair_coeff * * ${runnerCutoff}        # set up pair style coefficients

variable nnDir1    string "../nnp1/"
variable nnDir2    string "../nnp2/"
pair_style hybrid/overlay nnp dir ${nnDir1} showew no showewsum 10 resetew yes maxew 100000 cflength 1.8897261328 cfenergy 0.0015936014 nnp dir ${nnDir2} showew no showewsum 20 resetew yes maxew 100000 cflength 1.8897261328 cfenergy 0.0015936014
pair_coeff * * nnp 1 ${runnerCutoff}        # set up pair style coefficients
pair_coeff * * nnp 2 ${runnerCutoff}        # set up pair style coefficients
