# Response Envelope Analysis (REA)

This analytical tool enables the determination of drug combination effects (i.e., synergy, additivity, or antagonism) at local, regional, and global levels without knowing the inhibition mechanism of drugs. 

## Tutorial

The example file is "REA_example_run.m", which performs the analysis over representative examples from oncology and infectious disease settings. It calls the main routine function "REA_package.m". The other ".m" files are subroutines that are called by the main routine to estimate the Hill parameters and plot response envelope. The ".csv" files are the example data sets for the analysis. Please use the format for other data that need to analyzed using REA. 

## References
1. Du, D., et al, (2017), submitted
2. Griner, L.A.M., et al, High-throughput combinatorial screening identifies drugs that cooperate with ibrutinib to kill activated B-cell–like diffuse large B-cell lymphoma cells. Proc. Nat. Acad. Sci. (2014)
3. Borisy, A.A., et al, Systematic discovery of multicomponent therapeutics. Proc. Nat. Acad. Sci. (2003)
