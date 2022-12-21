-------------------------------------------------------------------------------
K-OPLS package for MATLAB
-------------------------------------------------------------------------------

Contents

1. Key features
2. Requirements
3. Installation
4. Getting started
5. License
6. Citing
7. Feedback 

-------------------------------------------------------------------------------
1. Key features
-------------------------------------------------------------------------------

The presented package provides an open-source, platform-independent
implementation of the Kernel-based Orthogonal Projections to Latent Structures
(K-OPLS) method; a kernel-based classification and regression method. In
relation to other kernel-based methods, K-OPLS offers unique properties
facilitating separate modeling of predictive variation and structured noise in
the feature space. While providing prediction results similar to other kernel-
based methods, K-OPLS features enhanced interpretational capabilities; allowing
detection of unanticipated systematic variation in the data such as
instrumental drift, batch variability or unexpected biological variation.

The package includes the following functionality:

(1.1)	Estimation (training) of K-OPLS models.

(1.2)	Prediction of new data using the estimated K-OPLS model in step (1.1). 

(1.3)	Cross-validation functionality to estimate the generalization error of
		a K-OPLS model. This is intended to guide the selection of the number
		of Y-predictive components A and the number of Y-orthogonal components
		Ao. The supported implementations are: 
		* n-fold cross-validation.
		* Monte Carlo Cross-Validation (MCCV)
		* Monte Carlo Class-balanced Cross-Validation (for discriminant 
		analysis cases).
	   
(1.4)	Kernel functions, including the polynomial and Gaussian kernel
		functions.
 
(1.5) 	Model statistics: 
		* The explained variation of X (R2X). 
		* The explained variation of Y (R2Y).
		* Prediction statistics over cross-validation for regression tasks
		(Q2Y, which is inversely proportional to the generalization error).
		* Prediction statistics over cross-validation for classification tasks 
		(sensitivity and specificity measures). 
	  
(1.6)	Plot functions for visualization: 
		* Scatter plot matrices for model score components.
		* Model statistics and diagnostics plots.

(1.7)   Automated optimization of the kernel parameter with either:
                * Gridsearch based on user-defined settings
                * Simulated annealing scheme


The K-OPLS package for MATLAB is freely available for download at
http://kopls.sourceforge.net/

-------------------------------------------------------------------------------
2. Requirements
-------------------------------------------------------------------------------

A functional installation of MATLAB 7.0 or later is required.

To get full functionality for the 'koplsPlotScores()' function, the
MATLAB statistics toolbox is required (see help on this function for further
details).


-------------------------------------------------------------------------------
3. Installation
-------------------------------------------------------------------------------

Unzip the MATLAB source files (*.m) into any directory of choice. Inside
MATLAB, set the path to include the extracted files.

-------------------------------------------------------------------------------
4. Getting started
-------------------------------------------------------------------------------

The package includes a demonstration of the functionality available in the
package, which is intended to be a starting point for future use. The
demonstration is based on a simulated data set, represented by 1000 spectral
variables from two different classes and is available in a supplied workspace.
The demonstration essentially consists two main steps.

The first step is to demonstrate how K-OPLS handles the model evaluation
(using cross-validation), model building and subsequent classification of
external data from a non-linear data set.

The second step is to demonstrate how K-OPLS works in the presence of response-
independent (Y-orthogonal) variation, using the same data set but with a strong
systematic class-specific disturbance added.

To run the demonstration, perform the following steps:
(1) Make sure that all included 'kopls*.m' files are in the current path.
(2) Run the demo by typing 'koplsDemo'.

-------------------------------------------------------------------------------
5. License
-------------------------------------------------------------------------------

The K-OPLS package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License version 2
as published by the Free Software Foundation.

The K-OPLS package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public
License version 2 along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.

-------------------------------------------------------------------------------
6. Citing
-------------------------------------------------------------------------------

When using this software in publications, please cite:

Rantalainen M, Bylesjö M, Cloarec O, Nicholson JK, Holmes E and Trygg J.
 Kernel-based orthogonal projections to latent structures (K-OPLS),
 J Chemometrics 2007; 21:376-385. doi:10.1002/cem.1071.

-------------------------------------------------------------------------------
7. Feedback
-------------------------------------------------------------------------------

Please send comments to Max Bylesjö <max_bylesjo@users.sf.net> or
Mattias Rantalainen <m_rantalainen@users.sf.net>
