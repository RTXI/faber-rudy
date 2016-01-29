###Faber-Rudy 2000

**Requirements:** None  
**Limitations:** None  

![Faber-Rudy GUI](faber-rudy.png)

<!--start-->

Faber-Rudy guinea pig cardiomyocyte model. 

**Reference:** <a href="http://dx.doi.org/10.1016/S0006-3495(00)76783-X">Faber GM, Rudy Y. "Action potential and contractility changes in Na+i overloaded cardiac myocytes: a simulation study". Biophys J. 2000 May;78(5):2392-404. Livshitz LM, Rudy Y. "Uniqueness and stability of action potential models during rest, pacing, and conduction using problem-solving environment". Biophys J. 2009 Sep 2;97(5):1265-76.</a>  

<!--end-->

####Input Channels
1. input(0) - Iapp : applied current  (A)

####Output Channels
1. output(0) - Vm : membrane voltage (V)

####Parameters
1. rate - rate of integration (Hz)
2. gNa - sodium conductance

####States
1. INa - sodium current
