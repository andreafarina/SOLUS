%% FIT4PARAM EXAMPLE

In this example we show how to reconstruct the coefficients of a inclusion from data simulated with TOAST++.

The flag TOAST2DOT in Init_DOT is set to 1 to allow conversion from asimulated data to a the format used in the software.

The extracted data undergo a convolution with experimental IRF and addition of Poisson noise. They are then reconstructed using a FEM implementation of the TD-DOT equation as forward problem.

The mask used from simulation and reconstruction is 'benign_1'.

To run this example, cd matlab to the example folder and run the run_DOT script.
