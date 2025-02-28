 #GXGD Distribution

The \textbf{GXGD (Generalized XGamma-Geometric Distribution)} is defined by the following probability density function (PDF):

\[
f_X(x; \theta, a, b, w) = \frac{w}{1 + w} \cdot \theta e^{-\theta x} + \frac{1}{1 + w} \cdot \frac{b^a}{\Gamma(a)} x^{a-1} e^{-b x}, \quad x > 0, \, \theta > 0, \, a > 0, \, b > 0, \, w > 0,
\]

where:
\begin{itemize}
    \item $\frac{w}{1 + w}$ and $\frac{1}{1 + w}$ are the mixing weights for the exponential and gamma components, respectively.
    \item $\theta e^{-\theta x}$ is the PDF of the exponential distribution.
    \item $\frac{b^a}{\Gamma(a)} x^{a-1} e^{-b x}$ is the PDF of the gamma distribution.
    \item $\Gamma(a)$ is the complete gamma function.
\end{itemize}

The cumulative distribution function (CDF) of the GXGD distribution is:

\[
F_X(x; \theta, a, b, w) = \frac{w}{1 + w} \left(1 - e^{-\theta x}\right) + \frac{1}{1 + w} \cdot \frac{\gamma(a, b x)}{\Gamma(a)}, \quad x > 0,
\]

where $\gamma(a, b x)$ is the lower incomplete gamma function.
 
