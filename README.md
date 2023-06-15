# Designing two-dimensional limit-cycle oscillators with prescribed trajectories and phase-response characteristics

__\<Description\>__ <br />
The code designs the vector field of two-dimensional limit-cycle oscillators with prescribed trajectories and phase-response characteristics.

You will get the figure like this:
<!--
<div align="center">
    <img src="figs/sl_phase.png">
</div>
<div style="text-align: center;">
    (a) True, (b) Estimation
</div> <br />
-->

This code is described in the paper
<!-- <a href="XXX" target="_blank">　-->
Designing two-dimensional limit-cycle oscillators with prescribed trajectories and phase-response characteristics</a> 
accepted by *XXX*. <br />

__\<Usage\>__ <br />
Please run 
- main.m for existing oscillators. 
- .m for the oscillator exhibiting multistable entrainment. 
- .m for the star-shaped oscillator. <br />

Files:
- floquet_2D.m: Calculation of PSF
- funcs.m: Functions for calculation (This contains "phase_numerical" function that obtains phase function numerically)
- main.m: File to run
- phase_response.m: Calculation of PRF
- polynomial_interpolation.m: Calculation of derivative of noisy data
- stuart_landau.m: Stuart-Landau oscillator
- FitzHugh-Nagumo.m: FitzHugh-Nagumo oscillator
- CIMA.m: Chlorite-iodide-malonic acid oscillator
- utils.m: utilities


__\<Environments\>__ <br />
MATLAB R2022a
- Optimization Toolbox
- Statistics and Machine Learning Toolbox

Python 3.9.6 (maybe work in other version)
- requirements.txt is added