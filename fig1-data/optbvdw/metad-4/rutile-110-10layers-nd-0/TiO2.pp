variable runnerCutoff    equal  6.352  # largest symmetry function cutoff (Angstrom)
variable runnerDir    string "/lustre1/u/yuechen/zengzz/2023/TiO2-water-defects/MLPs/optb88vdw/nnp-4"

mass 3 47.86
mass 2 16.00
mass 1 1.00

pair_style hdnnp ${runnerCutoff}  dir ${runnerDir} showew no showewsum 10 resetew yes maxew 500 cflength 1.8897261328 cfenergy 0.0015936014
pair_coeff * * H O Ti       # set up pair style coefficients
