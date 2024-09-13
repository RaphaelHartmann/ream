# ream
Calculate the probability density functions (PDFs) for two threshold evidence accumulation models (EAMs). These are defined using the following Stochastic Differential Equation (SDE), 

$`dx(t) = v(x(t),t)*dt+D(x(t),t)*dW,`$ 

where $`x(t)`$ is the accumulated evidence at time $`t`$, $`v(x(t),t)`$ is the drift rate, $`D(x(t),t)`$ is the noise scale, and $`W`$ is the standard Wiener process. The boundary conditions of this process are the upper and lower decision thresholds, represented by $`b_u(t)`$ and $`b_l(t)`$, respectively. Upper threshold $`b_u(t) > 0`$, while lower threshold $`b_l(t) < 0`$. The initial condition of this process $`x(0) = z`$ where $`b_l(t) < z < b_u(t)`$. We represent this as the relative start point $`w = z/(b_u(0)-b_l(0))`$, defined as a ratio of the initial threshold location. This package generates the PDF using the same approach as the Python package it is based upon, PyBEAM by Murrow and Holmes (2023) <doi:10.3758/s13428-023-02162-w>. First, it converts the SDE model into the forwards Fokker-Planck equation 

$`dp(x,t)/dt = d(v(x,t)*p(x,t))/dt-0.5*d^2(D(x,t)^2*p(x,t))/dx^2,`$

then solves this equation using the Crank-Nicolson method to determine $`p(x,t)`$. Finally, it calculates the flux at the decision thresholds, $`f_i(t) = 0.5*d(D(x,t)^2*p(x,t))/dx`$ evaluated at $`x = b_i(t)`$, where $`i`$ is the relevant decision threshold, either upper ($`i = u`$) or lower ($`i = l`$). The flux at each thresholds $`f_i(t)`$ is the PDF for each threshold, specifically its PDF. We discuss further details of this approach in this package and PyBEAM publications. Additionally, one can calculate the cumulative distribution functions of and sampling from the EAMs.
