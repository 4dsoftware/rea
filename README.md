# Response Envelope Analysis (REA)

This analytical tool enables the determination of drug combination effects (i.e., synergy, additivity, or antagonism) at local, regional, and global levels without knowing the inhibition mechanism of drugs. 

## Instruction

The example file is "REA_example_run.m", which performs the analysis over representative examples from oncology and infectious disease settings. It calls the main routine function "REA_package.m". The other ".m" files are subroutines that are called by the main routine to estimate the Hill parameters and plot response envelope. The ".csv" files are the example data sets for the analysis. Please use the format for other data that need to analyzed using REA. 

## Tutorial

To use the REA package for other applications, one can simply call the function REA_package.m. The m-function requires four input variables, data, trim, ft, and lw. "data" is the data array obtained from experiments or simulations. The format is given in the example and should be used. "trim" is the margin of the plot beyond exisiting measurements. "ft" and "lw" are font size and line width, respectively. The last three are rather aesthetical. Other optional input variables include drug names and "custom_label" which specifies whether custom axis labels are used. If not customized, logarithm of axis labels will be shown. 

## References
1. Du, D., et al, (2017), submitted
2. Griner, L.A.M., et al, High-throughput combinatorial screening identifies drugs that cooperate with ibrutinib to kill activated B-cell–like diffuse large B-cell lymphoma cells. Proc. Nat. Acad. Sci. (2014)
3. Borisy, A.A., et al, Systematic discovery of multicomponent therapeutics. Proc. Nat. Acad. Sci. (2003)
