Proposed outline:

1. Introduction
2. Parametric models
    * Review the general theory of parametric inversion (like in text books)
    * This will serve to establish the naming conventions and mathematical
      notation
    * Show only the fully non-linear case
    * Include Tikhonov regularization
3. Classical theory of implicit models
    * Review the derivations in the Geodesy: the concepts book
    * Show equations for the simplest case (only minimize norm of data step)
    * Show equations for case minimizing the parameter step as well
    * Mention that minimizing delta p is like using Levemberg-Marquardt, not
      regularization (compare with Tikhonov 0 expression deduced above)
4. A new look at implicit models for geophysics
    * Highlight some of the strange points in the classical equations
        * Minimizing the step in the predicted data, not the residuals
        * No observed data in the equations
        * Objective function is formulated after linearization
        * Need proper regularization
    * Formulate the problem using all available tools:
        * Norm of residuals
        * Tikhonov regularization
        * Norm of data step
        * Norm of parameter step
5. Example from geophysics
    * The best candidate is Love wave dispersion with a layer + half space
      model
    * Other options: Thermobarometry, Euler deconv
6. Conclusions
