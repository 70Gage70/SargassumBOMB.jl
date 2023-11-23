# The eBOM equations

Refer to [here](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/minimal-maxeyriley-model-for-the-drift-of-sargassum-rafts/5236FED7891D4309789EF0696907655F).

Particles have radius $a > 0$ and water-to-particle density $\delta \geq 1$. We define
$$
\begin{align}
\psi(\delta) &= \sqrt[3]{i \sqrt{1-\left(\frac{2}{\delta }-1\right)^2}+\frac{2}{\delta }-1}, \\ 
\Phi(\delta) &= \frac{1}{2} \left(i \sqrt{3}\right) \left(\frac{1}{\psi (\delta )}-\psi (\delta )\right)-\frac{1}{2} \left(\psi (\delta )+\frac{1}{\psi (\delta )}\right)+1, \\ 
\Psi(\delta) &=  \pi^{-1} \cos ^{-1}(1-\Phi (\delta )) -  \pi^{-1} (1-\Phi (\delta )) \sqrt{1-(1-\Phi (\delta ))^2}.
\end{align}
$$
Let $\rho$ be the (constant) density of the water and $\mu$ be the water dynamic viscosity. Let $\gamma$ be the air-to-water viscosity ratio. Then,
$$
\begin{align}
\alpha(\delta) &= \frac{\gamma \Psi}{1 - (1 - \gamma) \Psi}, \\
\tau(\delta) &= \frac{1 - \tfrac{1}{6} \Phi}{(1 - (1 - \gamma) \Psi) \delta^4} \frac{a^2 \rho}{3 \mu} , \\ 
R(\delta) &= \frac{1 - \tfrac{1}{2} \Phi}{1 - \tfrac{1}{6} \Phi} .
\end{align}
$$
Let $\Omega$ be the angular velocity of the Earth in radians and then $\vartheta_0$ be a reference latitude. Then,
$$
\begin{equation}
f = 2 \Omega \sin \vartheta_0 .
\end{equation}
$$
Let the planar $(x, y)$ velocity of the water be $v(x, y, t)$, let the planar velocity of the air be $v_{a}(x, y, t)$ and define
$$
\begin{equation}
u = (1 - \alpha) v + \alpha v_a.
\end{equation}
$$
Let $\omega(x, y, t)$ be the vorticity of the water, then
$$
\begin{equation}
\omega = \frac{\partial v_y}{\partial x} - \frac{\partial v_x}{\partial y}.
\end{equation}
$$

Let $\text{D}v_i/\text{D}t = (v, 1) \cdot (\nabla_{x, y, t} v_i) = \partial_t v_i + v \cdot (\nabla_{x, y} v_i) $.
$$
\begin{align}
\dot{x} &= u_x + \tau\left( R \frac{\text{D} v_x}{\text{D}t} - R\left( f + \frac{1}{3}\omega \right)v_y -  \frac{\text{D} u_x}{\text{D}t} + \left( f + \frac{1}{3}R \omega \right) u_y  + F_x\right) \\ 
\dot{y} &= u_y + \tau\left( R \frac{\text{D} v_y}{\text{D}t} + R\left( f + \frac{1}{3}\omega \right)v_x -  \frac{\text{D} u_y}{\text{D}t} - \left( f + \frac{1}{3}R \omega \right) u_x  + F_y \right)
\end{align}
$$