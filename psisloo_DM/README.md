# psisloo
PSISLOO by Aki Vehtari et al. with small modifications to avoid rare errors.
 
*Pareto-smoothed importance sampling leave-one-out cross-validation* (PSISLOO) is a model comparison measure developed by Aki Vehtari et al.: https://github.com/avehtari/PSIS. 
This is a slight modification of their Matlab code. 

I made a very small change in psislw.m: catching some NaNs that would otherwise result in an erroneous NaN PSISLOO outcome. 
Another small change was made in gpdfitnew.m: avoiding division by zero.
