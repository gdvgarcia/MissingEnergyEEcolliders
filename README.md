# Electron-positron Colliders Event Generator for n(iDM) Package

## What is it?

This package allows the user to generate a variety of high energy \
events at electron-positron colliders in the context of the Standard \
Model supplemen-
\n ted by not-so-inelastic Dark Matter model (DM) with a dark photon \
mediator. DM notation will refer to the lighest dark fermionic state \
\\ chi- (\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)), e.g. 
\nmDM = \!\(\*SubscriptBox[\(m\), \(\[Chi]\)]\) = mass of chi-. Zp \
(Z') will refer to the dark photon.

This package is a Mathematica-based code allowing to calculate constraints/sensitivities from electron-positron colliders for (n)iDM models. The method of the calculation is described in [2305.13383](https://arxiv.org/abs/2305.13383). 


## How to launch

Download the full directory and check the example code included.

### Dependencies

To run our package, two tools have to be installed: Wolfram Mathematica (tested on version 13.1) and a C compiler. 


## (Currently) implemented

For the moment, only the usual niDM model with a dark photon mediator is implemented (implementing other models would require substantial work). 

Also, we have only focused on missing energies so far. However, for other kinds of searches, one can simply compute the total number of events ab inition using the total cross-section and the integrated luminosity of the detector, then, normalize the number of events simulated.

Last, we currently only focused on BaBar and Belle II. Although, considering other e+e- colliders should not be a difficult task.
