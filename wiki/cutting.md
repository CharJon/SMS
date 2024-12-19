# SCIPs cut-pool management / cut selection

SCIPs cut selection has two phases, a scoring and a selection phase

## Cut scoring

Every cut gets scored based on:

- Efficacy, i.e., distance of the hyperplane to the LP sol: $`e_r`$
- Orthogonality with respect to the other cuts: $`o_r`$
- Parallelism with respect to the objective function: $`p_r`$

The final formula for a cuts score is just: $`e_r + w_o * o_r + w_p * p_r`$.
- Default values in SCIP for $`w_o`$ named "ORTHOFAC" is 1.
- Default values in SCIP for $`w_p`$ named "OBJPARALFAC" is 0.0001.

## Cut selection

After every cut got scored, the cutpool gets sorted based on the cuts scores.
Let $`c_1`$ be the maximum score. SCIP assigns every cut into one of three categories: good, decent and bad.
A cut with a score >= "goodscorefac" * $`c_1`$ is considered good (default for "goodscorefac" is 0.9).
A cut with a score <= "badscorefac" * $`c_1`$ is considered bad (badscorefac is 0 by default, which means this is turned off in default settings).
Every other cut is considered to be decent.

Cuts in the "bad" category get disregarded completely.
To finally select the cuts SCIP loops over all cuts in non-increasing order of their score
and picks a cut if it is at most $`max(0.5,n)`$ parallel to _all_ previously selected cuts if it is a "good" cut
and at most $`1-n`$ parallel to _all_ previously selected cuts if it is a "decent" cut.
$`n`$ named "maxparall" in SCIP has a default value of 0.9.
