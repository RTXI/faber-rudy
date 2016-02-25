###Faber-Rudy 2000

**Requirements:** None  
**Limitations:** None  

![Faber-Rudy GUI](faber-rudy.png)

<!--start-->

The Faber-Rudy guinea pig cardiomyocyte model. Scaling parameters in the UI allow you to tweak model parameters within RTXI without editing the source code and recompiling. 

**Reference:** 

 * <a href="http://dx.doi.org/10.1016/S0006-3495(00)76783-X">Faber GM, Rudy Y. "Action potential and contractility changes in Na+i overloaded cardiac myocytes: a simulation study". Biophys J. 2000 May;78(5):2392-404.</a>
 * <a href="http://dx.doi.org/10.1016/j.bpj.2009.05.062">Livshitz LM, Rudy Y. "Uniqueness and stability of action potential models during rest, pacing, and conduction using problem-solving environment". Biophys J. 2009 Sep 2;97(5):1265-76.</a>  

<!--end-->

####Input Channels
1. input(0) - Iapp : applied current  (A)

####Output Channels
1. output(0) - Vm : membrane voltage (V)

####Parameters
1. rate - rate of integration (Hz)
2. gNa - sodium conductance
3. Scale GNa - scaling factor for sodium conductance
4. Scale GCaL - scaling factor for L-type calcium conductance
5. Scale GCaT - scaling factor for T-type calcium conductance
6. Scale GK1 - scaling factor for IK1 current
7. Scale GKr - scaling factor for IKr current
8. Scale GKs - scaling factor for IKs current
9. Scale GKp - scaling factor for IKp current
10. Scale IpCa - scaling factor for IpCa
11. Scale JSERCA - scaling factor for JSERCA

####States
1. INa - sodium current
