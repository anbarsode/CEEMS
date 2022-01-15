# CEEMS  
**C**oupled **E**quivalent **E**qual **M**ass **S**ystems for classical central force problems.  
A generalization of the classic equivalent one body problem.  

![CEEMS_random](./CEEMS_random.png "Demonstration with random inputs")  

This program illustrates the CEEMS formalism via numerical simulations. The Lagrangian for motion of two particles of mass under a mutual central potential is given by
<p align="center"> <img src="https://render.githubusercontent.com/render/math?math=\mathcal{L}=\dfrac{1}{2}m_1\dot{r}^2_1%2B\dfrac{1}{2}m_2\dot{r}^2_2-U(\lvert \mathbf{r}_1 - \mathbf{r}_1 \rvert)"></p>  

We would like to obtain an **E**quivalent Lagrangian of the form  
<p align="center"> <img src="https://render.githubusercontent.com/render/math?math=\mathcal{L}=\dfrac{1}{2}A\dot{r}^2_a%2B\dfrac{1}{2}A\dot{r}^2_b-U(\lvert \mathbf{r}_a - \mathbf{r}_b) \rvert/\gamma)"></p>  
with  
<p align="center"> <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r}_a=\alpha_1\mathbf{r}_1%2B\alpha_2\mathbf{r}_2"></p>  
<p align="center"> <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r}_b=\beta_1\mathbf{r}_1%2B\beta_2\mathbf{r}_2"></p>  
where <img src="https://render.githubusercontent.com/render/math?math=A"> is some **M**ass, **E**qual for both the particles, and it, along with <img src="https://render.githubusercontent.com/render/math?math=\alpha_1, \alpha_2, \beta_1, \beta_2, \gamma_1, \gamma_2"> are constants that can only depend on the two masses. Thus the relationship between the two systems is linear, the potential is still central, and the two masses are **C**oupled (unlike the classic equivalent one body problem).  

The interpretation of this new system is as follows: Instead of obtaining equations of motion using Euler-Lagrange equations with the original Lagrangian, under a potential <img src="https://render.githubusercontent.com/render/math?math=U(\lvert \mathbf{r}_1 - \mathbf{r}_1 \rvert)">, with initial conditions <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r}_1(t_0),\mathbf{\dot{r}}_1(t_0),\mathbf{r}_2(t_0),\mathbf{\dot{r}}_2(t_0)">; we obtain equations of motion using the equivalent Lagrangian, under the same functional form of the potential, but with different coefficients (even those are identical in two of the models) and with different initial conditions <img src="https://render.githubusercontent.com/render/math?math=\mathbf{r}_a(t_0)=\alpha_1\mathbf{r}_1(t_0)+\alpha_2\mathbf{r}_2(t_0), \mathbf{\dot{r}}_a(t_0)=\alpha_1\mathbf{\dot{r}}_1(t_0)+\alpha_2\mathbf{\dot{r}}_2(t_0), \mathbf{r}_b(t_0)=\beta_1\mathbf{r}_1(t_0)+\beta_2\mathbf{r}_2(t_0), \mathbf{\dot{r}}_b(t_0)=\beta_1\mathbf{\dot{r}}_1(t_0)+\beta_2\mathbf{\dot{r}}_2(t_0)">.  

By equating the two Lagrangians by vectorial terms, a system of algebraic equations can be obtained and solved to determine the transformation parameters. In the three models coded in `CEEMS.py`, one corresponds to the case where a value for the equal mass is randomly chosen and the equivalent system evolves under the action of a force with scaled coefficients. The other two models correspond to the cases where the force coefficients are kept the same by a suitable choice of the equivalent equal mass. For systems where the force coefficients do not depend on mass (anything other than gravity basically), the choice is <img src="https://render.githubusercontent.com/render/math?math=A=2\mu">. For the case of Newtonian gravity, the choice is made so as to keep the value of G same in the equivalent "Universe". This leads to <img src="https://render.githubusercontent.com/render/math?math=A=2^{1/5}\mathcal{M}"> (proportional to the Chirp mass of the binary system). Incredibly, this choice directly causes the leading order gravitational wave emission by these systems to be identical in amplitude, frequency and energy flux for **any** mass ratio. This is also one of the reasons why the CEEMS formalism doens't hold in General Relativity.  

Do let me know if you find any aplications for this.

Edit model parameters in `CEEMS.py`.  

To run the code,  
```
$ py CEEMS.py
```
