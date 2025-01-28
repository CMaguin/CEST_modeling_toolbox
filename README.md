# CEST_modeling_toolbox

## Summary

This is a matlab-based toolbox for CEST modeling and quantitative fitting.

This toolbox is associated with my PhD work, entitled "Modeling of CEST-MRI signal for the quantification of brain metabolism", defended on the 10th of October 2024. 
For documentation, you can have a look at the PhD manuscript (unfortunately only available in French) which can be found here : 
[https://theses.hal.science/tel-04810545](https://theses.hal.science/tel-04810545)

As a heads up, the codes found here are to some extent inspired by the codes of [cest-sources of M.Zaiss lab](https://github.com/cest-sources). You can also find there a lot of useful codes and documentation to understand and improve CEST modeling. 

## Examples of use

Some examples to use this code can be found in "Examples" folder. A toy data set of Glutamate phantoms data is also provided. Details about the data are included in the read-me of the "ToyDataSet" folder.

You can just run the scripts after adding all codes to path (or run the console command addpath(genpath(pwd)) ).

### Example 1 : Numerical versus pseudo-analytical simulation

With the "Comparison_NumAnaSimulations.m" code, you can find an example of how to quickly compute CEST simulations with your typical experimental setup with either numerical solution or pseudo-analytical solution (see [Zaiss & Bachert 2013](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/nbm.2887) ) and assess the time it takes to compute both. 

In short, numerical simulations are more rigorous, but take typically 10-100 more times to run, and can be especially long with complex saturation pulse shape design. Pseudo-analytical simulations are of course faster to compute, but can only handle some extent of complexity, and can deviate a bit from "ground-truth" in some limit cases (see my manuscript, Chapter IV for more details). Depending on your CEST design and numerical ressources, you might choose to compromise for one or the other solution. 

![Comparison of numerical and analytical simulations](/Examples/illustrating-figures/Fig1.png)

### Example 2 : some insightful simulations

Some examples of CEST simulations which can provide some insights into understanding the complexity of the CEST phenomenon. These are mainly associated with Chapter III of my thesis, dealing with the question of design of a CEST experiment. 
Here is a non-exhaustive list :

+ **Simulating how an experimental parameter affects the CEST signal** : For instance, you can look at how the B<sub>0</sub> field strength affects the CEST signal ("Simulation_B0influence.m" and "Simulation_B0influence_exogeneous.m"). You can understand how increasing B<sub>0</sub> strength migth be useful for CEST agents close to water resonance (endogeneous agents for example), but might be overkill for exogeneous CEST agents, far from water resonancy.

+ **Optimize saturation parameters to maximize the CEST contrast of your agent of interest** : With the script "Sensitivity_map_example_Glu_1000Hz.m", you can see how you can generate what I call a "sensivity-map". The intensity of the CEST contrast is plotted in a color-map as a function of B<sub>1</sub> and t<sub>sat</sub> (CW case here, but you can also try pulse). This can help you optimize your saturation module parameters. The optimized parameters are dependent of your B<sub>0</sub> field, but are also heavily dependent on the assumption of the exchange rate value k<sub>ex</sub>, which is often still an uncertain value in the CEST research community.

![Sensitivity map at 11.7T for glutamate with 1000 Hz exchange rate](/Examples/illustrating-figures/Fig2.png)
  
+ **Simulating how the CEST signal would look like in vivo**, for instance at 11.7T in my case (script "Simulation_Zspectrum_invivo117T.m") : you can test out different combinations of included CEST agents and vary experimental parameters. This can help you assess how much your CEST signal of interest could become contaminated by the presence of other endogeneous CEST agents in vivo. This might help you design how you want to sample saturation offsets during your real CEST experiments, depending on what you are interested in measuring/removing. Careful though, simulations might not always correspond perfectly to reality because of the great number of parameters involved and their associated uncertainties. As best practice, it is preferable to confirm what your Z-spectrum looks like in vivo with a densely sampled range of saturation frequencies.

+ **Simulating more quantitatively how contaminated your CEST contrast is** : with the script "Simulating_Contributions_to_MTR_based_on_KhlebnikovEtAl.m", you assess more accurately how much your agent of interest contribute to the actual CEST signal expected in vivo. The method and chose parameters were taken from [Khlebnikov et al](https://www.nature.com/articles/s41598-018-37295-y). Here is an illustration below.
Although the adding up of individual CEST contrasts can pretty accurately predict the proportions of contributions of CEST agents, note that the actual resulting CEST contrast is not as high as the addition of individual CEST contast (this is because the CEST phenomenon is non linear).

![Expected contributions of metabolites to in vivo MTR](/Examples/illustrating-figures/Fig3.png)

### Example 3 : Fitting toolbox

With the script "Fitting_exchange_rate_GluCEST_ToyDataSet.m", I propose a fitting function to estimate exchange rate values from multi-B<sub>1</sub> experimental CEST data (namely the glutamate toy data set). You can look at how the buffer affects the exchange rate value. You can also expect that the actual exchange rate value would be vastely different in vivo than from in vitro. 

![Glutamate phantoms at 37C in different buffers](/Examples/illustrating-figures/Fig4.png)  ![Fitting multiB1 data](/Examples/illustrating-figures/Fig5.png)


