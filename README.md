# Response Envelope Analysis (REA)

This analytical tool enables the determination of drug combination effects (i.e., synergy, additivity, or antagonism) at local, regional, and global levels without knowing the inhibition mechanism of drugs. It was developed under MATLAB computing environment, and all the scripts are ".m" files. Some useful variants for users with different background are down below: 

1. REA with graphic interface (executable, no scripting is needed), https://github.com/4dsoftware/grea. 
2. Three-drug combination, https://github.com/4dsoftware/rea43.
3. REA on R, https://github.com/4dsoftware/reaR. __NEW__!

![cover](https://user-images.githubusercontent.com/15344717/32510711-182d3990-c3b7-11e7-9fbc-b1d796fc3706.jpg)

## Instruction

The example file is "REA_example.m", which performs the analysis over representative examples from oncology and infectious disease settings. It calls the main routine function "REA_package.m". The other "package_*.m" files are subroutines that are called by the main routine to estimate the Hill parameters and plot response envelope. The ".csv" files are the example data sets for the analysis. Please use the format for other data that need to be analyzed using REA. Please also note that a zero dose and at least 3-by-3 dose matrix are required to run the application. 

## Tutorial

To use the REA package for other applications, one can simply call the function "REA_package.m". The m-function requires four input variables, data, trim, ft, and lw. "data" is the data array obtained from experiments or simulations. The format is given in the example and should always be used. "trim" is the margin of the plot beyond exisiting measurements. "ft" and "lw" are font size and line width, respectively. The last three are rather aesthetical. Other optional input variables include drug names and "custom_label" which specifies whether custom axis labels are used. If not customized, logarithm of axis labels will be shown. 

The output of the m-function "REA_package.m" includes the synergy index and antagonism index and the pipeline used to obtain these two indexes. The synergy index and antagonism index reflect regional combination effects, and the difference between these two suggests global effects.

Please see "REA_Tutorial.pdf" for more details.

## References
1. Du, D., et al, Response Envelope Analysis for Quantitative Evaluation of Drug Combinations, Bioinformatics, (2019)
2. Griner, L.A.M., et al, High-throughput combinatorial screening identifies drugs that cooperate with ibrutinib to kill activated B-cell–like diffuse large B-cell lymphoma cells. Proc. Nat. Acad. Sci. (2014)
3. Borisy, A.A., et al, Systematic discovery of multicomponent therapeutics. Proc. Nat. Acad. Sci. (2003)

## Contact
Dr D. Du, dudthu06@gmail.com
