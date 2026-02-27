# Tukey Additive Spatial Model

## 📝 Description
Implementation of an MCMC algorithm for a spatial model with Tukey transformation for heavy-tailed data.

## ⚙️ Important Note
**Documentation for function parameters is provided within the function itself.**
**These objects must be set before running the function.**

## 📋 Model Parameters
The following parameters need to be set in the code before execution:
- `phi`: Spatial decay parameter
- `sig`: Scale parameter for V component  
- `tau2`: Error variance
- `h`: Tail heaviness parameter
- `alpha`: Skewness parameter
- `beta`: Regression coefficient
- `U`, `V`: Initial values for latent spatial processes
- tun: Tunning parameter
- n.n:  Number of neighbors

## 🚀 Usage
```r
# Set your parameters first, then run:
Tukey(simulate = 10000, burnin = 2000, Break = 10, tun = 0.5)
