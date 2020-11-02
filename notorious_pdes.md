# Notorious PDEs

These PDEs are generally regarded by the community (or by me in some cases) as being quite difficult to nearly impossible to solve analytically. Simulating these numerically is my long-term goal.

- Primitive Equations
    - ![$\begin{aligned}\frac{dv}{dt}&=-(1/\rho)\nabla p-g(r/r)+f_{r}\\ g&=g_e\\ \frac{dv}{dt}&=-(1/\rho)\nabla p-g(r/r)+(1/\rho )\left[\nabla\cdot (\mu\nabla v)+\nabla (\lambda\nabla\cdot v)\right]\\ c_{v}\frac{dT}{dt}+p\frac{d\alpha}{dt}&=q+f\\ \frac{d\rho}{dt}+\rho\nabla\cdot v&=0\\ p&=nT\\ \end{aligned}$](svg/primitive_equations.svg)
- Einstein Field Equations
    - ![$R_{\mu\nu}-\frac{1}{2}Rg_{\mu\nu}+\Lambda g_{\mu\nu}=\frac{8\pi G}{c^4}T_{\mu\nu}$](svg/einstein_field_equations.svg)
- Navier-Stokes Equations
    - Convective Form
        - ![$\rho\frac{D\mathbf{u}}{Dt}=\rho\left(\frac{\partial\mathbf{u}}{\partial t}+\mathbf{u}\cdot\nabla\mathbf{u}\right)=-\nabla\bar{p}+\mu\,\nabla^{2}\mathbf{u}+\frac{1}{3}\mu\,\nabla(\nabla\cdot\mathbf{u})+\rho\mathbf{g}$](svg/navier_stokes_equations_convective.svg)
    - Conservation Form
        - ![$\frac{\partial}{\partial t}(\rho\,\mathbf{u})+\nabla\cdot(\rho\,\mathbf{u}\otimes\mathbf{u})=-\nabla\bar{p}+\mu\,\nabla^{2}\mathbf{u}+\frac{1}{3}\mu\,\nabla(\nabla\cdot\mathbf{u})+\rho\mathbf{g}$](svg/navier_stokes_equations_conservation.svg)
- Maxwell's Equations
    - ![$\begin{aligned}\nabla\cdot\mathbf{E}&=\frac{\rho}{\varepsilon_0}\\\nabla\cdot\mathbf{B}&=0
\\\nabla\times\mathbf{E}&=-\frac{\partial\mathbf{B}}{\partial t}
\\\nabla\times\mathbf{B}&=\mu_0\left(\mathbf{J}+\varepsilon_0\frac{\partial\mathbf{E}}{\partial t}\right)\\\end{aligned}$](svg/maxwells_equations.svg)
- Shallow Water Equations
    - Conservative Form
        - ![$\begin{aligned}\frac{\partial(\rho\eta)}{\partial t}+\frac{\partial(\rho\eta u)}{\partial x}+\frac{\partial(\rho\eta v)}{\partial y}&=0\\ \frac{\partial (\rho\eta u)}{\partial t}+\frac{\partial}{\partial x}\left(\rho\eta u^2+\frac{1}{2}\rho g\eta^2\right)+\frac{\partial(\rho\eta uv)}{\partial y}&=0\\ \frac{\partial(\rho\eta v)}{\partial t}+\frac{\partial(\rho\eta uv)}{\partial x}+\frac{\partial}{\partial y}\left(\rho\eta v^2+\frac{1}{2}\rho g\eta^2\right)&=0\end{aligned}$](svg/shallow_water_equations_conservative.svg)
    - Non-Conservative Form
        - ![$\begin{aligned}\frac{\partial h}{\partial t}+\frac{\partial}{\partial x}\Bigl((H+h)u\Bigr)+\frac{\partial}{\partial y}\Bigl((H+h)v\Bigr)&=0\\\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}+v\frac{\partial u}{\partial y}-fv&=-g\frac{\partial h}{\partial x}-bu+\nu\left(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\right)\\ \frac{\partial v}{\partial t}+u\frac{\partial v}{\partial x}+v\frac{\partial v}{\partial y}+fu&=-g\frac{\partial h}{\partial y}-bv+\nu\left(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\right)\\ \end{aligned}$](svg/shallow_water_equations_nonconservative.svg)
- Universal Differential Equation(s)
    - ?

### Honorable Mention
These might not be as ubiquitous or challenging to solve, but I still find them cool and want to simulate them. I might move them to the main list in the future...
- Lamb-Oseen Vortex
    - ![$\frac{\partial g}{\partial t}=\nu\left(\frac{\partial^2g}{\partial r^2}-\frac{1}{r}\frac{\partial g}{\partial r}\right)$](svg/lamb_oseen_vortex.svg)
- Thin-Film Equation
    - ![$\frac{\partial h}{\partial t}=-\frac{1}{3\mu}\nabla\cdot\left(h^3\,\nabla\left(\gamma\,\nabla^2h\right)\right)$](svg/thin_film_equation.svg)
- Vorticity Equation
    - ![$\begin{aligned}\frac{D\boldsymbol{\omega}}{Dt}&=\frac{\partial\boldsymbol{\omega}}{\partial t}+(\mathbf{u}\cdot\nabla)\boldsymbol{\omega}\\&=(\boldsymbol{\omega}\cdot\nabla)\boldsymbol{u}-\boldsymbol{\omega}(\nabla\cdot\boldsymbol{u})+\frac{1}{\rho^2}\nabla\rho\times\nabla p+\nabla\times\left(\frac{\nabla\cdot\tau}{\rho}\right)+\nabla\times\left(\frac{\boldsymbol{B}}{\rho}\right)\end{aligned}$](svg/vorticity_equation.svg)
