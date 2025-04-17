# FunctionalVAR

This repository implements the **Functional VAR (FVAR)** approach from Huang (2023). models the joint dynamics of macro aggregates and functional variables within the SVAR framework. It supports three identification schemes: (1) short-run restrictions, (2) internal instrumental variables (IV), and (3) external IV. Inference is performed via a moving-block bootstrap variant.

For theoretical details, refer to the [latest draft](https://jiaminghuang.net). For implementation questions, contact [jiaming.huang@barcelonagse.eu](mailto:jiaming.huang@barcelonagse.eu).

---

## Minimal Working Example

See `MWE.m` for minimal implementation of the simulation exercise in the paper. For your application with macro data (`Txn`) and functional data (`TxnGridFcn`), you can run the following:

```matlab
DATASET.data = data;                    % T x n matrix of macro variables
DATASET.gridFcn = gridOfYourFcn;        % 1 x ngridFcn vector of function grid points
DATASET.fcn = fcnData;                  % T x ngridFcn matrix of functional data
DATASET.varsName = colNamesOfYourData;  % Cell array of variable names
modelSpec.varsSel = varsToIncludeInVAR; % e.g., {'x', 'Func', 'y'}
FVAR = doFVAR(DATASET, modelSpec);      % Run estimation
```

The output `FVAR` contains:

- Estimated VAR coefficients
- Model specification
- Impulse response functions (IRs)
- Bootstrap confidence intervals

Then you can play with the model to data setup to your needs (see below for supported features).

---

## Specification

### 1. `DATASET` Struct

| Field      | Description                                                                  |
| ---------- | ---------------------------------------------------------------------------- |
| `data`     | `T x n` matrix of macro variables                                            |
| `varsName` | `1 x n` cell array of variable names                                         |
| `varsDiff` | Cell array of differenced variables (cumulative IRs computed). Default: `{}` |
| `fcn`      | `T x ngridFcn` matrix of functional data evaluations                         |
| `gridFcn`  | `1 x ngridFcn` vector of function domain grid points                         |
| `dates`    | `T x 1` vector of dates (`datetime` or numeric)                              |

### 2. `modelSpec` Struct

#### VAR Specification

| Parameter                | Description                                                                                                              |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------ |
| `varsSel`                | Cell array of variables in VAR order (Cholesky ordering if using identification='CHOL')<br>Example: `{'x', 'Func', 'y'}` |
| `p`                      | VAR lag order. Default: `1`                                                                                              |
| `iDET`                   | Deterministic components:<br>`1`=constant, `2`=linear trend, `3`=quadratic trend. Default: `1`                           |
| `dateStart`<br>`dateEnd` | Estimation sample time window                                                                                            |

#### FPCA Specification

| Parameter  | Description                                                                                                         |
| ---------- | ------------------------------------------------------------------------------------------------------------------- |
| `fdParobj` | `fdPar` object for functional smoothing (Default: Fourier basis):<br>```matlab create_fourier_basis([0,1], 21); ``` |
| `nFPCMax`  | Max FPCs to consider. Selected via cumulative variance. Default: `20`                                               |

#### SVAR Specification

| Parameter        | Description                                                           |
| ---------------- | --------------------------------------------------------------------- |
| `identification` | `"CHOL"` (Cholesky), `"InternalIV"`, or `"SVARIV"`. Default: `"CHOL"` |
| `varsShock`      | Shock variable(s). Default: First entry in `varsSel`                  |
| `varsUnitNorm`   | Variables for unit shock normalization (Internal IV only)             |
| `Instrument`     | Instrument variable (Required for `"SVARIV"`)                         |
| `nZlags`         | IV lag order (External IV only). Default: `0`                         |
| `nNWlags`        | Newey-West covariance lags. Default: `floor(1.3*sqrt(T))`             |
| `irhor`          | Maximum IRs horizon (total horizons = `irhor + 1`). Default: `36`     |

#### Bootstrap

| Parameter | Description                                  |
| --------- | -------------------------------------------- |
| `nBoot`   | Bootstrap replications. Default: `500`       |
| `blkSize` | Moving block length. Default: `16`           |
| `cLevel`  | Confidence level (%). Default: `95`          |
| `verbose` | Display bootstrap progress. Default: `false` |

---

## FPCA Implementation  

### FDA Package Requirement

I use the `fdaM` MATLAB package for Functional Principal Component Analysis (FPCA):  
**[Download fdaM](https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/)**.

### 1. Functional Data Smoothing

See Chapter 3 of Ramsay et al. (2009) for details on the basis functions. For example, the Fourier basis can be specified as follows:

```matlab
% Create Fourier basis
rangeval = [0, 1];       % Function domain
nbasis = 21;             % Number of basis functions
basisobj = create_fourier_basis(rangeval, nbasis);

% Convert raw data to functional objects
fdParobj = fdPar(basisobj);
fdobj = smooth_basis(argvals, f, fdParobj);  % f: nGrid x T matrix
```

### 2. FPCA

```matlab
% Perform FPCA
pcaOut = pca_fd(fdobj, modelSpec.nFPCMax, fdParobj);

% Extract components
mu = eval_fd(gridFcn, pcaOut.meanfd)';     % Mean function
bas = eval_fd(gridFcn, pcaOut.harmfd);    % Eigenfunctions (nGrid x nFPC)
fpcs = pcaOut.harmscr;                  % FPC scores (T x nFPC)

% Retain FPCs explaining ≥95% variance
nFPC = find(cumsum(pcaOut.varprop) >= 0.95, 1);
fhat = mu + fpcs(:,1:nFPC) * bas(:,1:nFPC)';  % Reconstructed functions
```

---

## SVAR Implementation

- **Short-run restrictions**: Recursive identification via Cholesky (Kilian & Lütkepohl, 2017).  
- **Internal IV**: Instrument ordered first in a recursive VAR (Plagborg‐Møller & Wolf, 2021).  
- **External IV**: Reduced-form VAR, plus a separate IV regression (Stock & Watson, 2018).  

**Bootstrap**:  

- Recursive SVAR: Brüggemann et al. (2016).  
- SVAR-IV: Jentsch & Lunsford (2022).  

---

## References

1. Brüggemann, R., Jentsch, C., & Trenkler, C. (2016). Inference in VARs with conditional heteroskedasticity of unknown form. *Journal of Econometrics*, *191*(1), 69–85.  
2. Huang, J. (2023). Functional VAR.  
3. Jentsch, C., & Lunsford, K. G. (2022). Asymptotically valid bootstrap inference for proxy SVARs. *Journal of Business & Economic Statistics*, *40*(4), 1876–1891.  
4. Kilian, L., & Lütkepohl, H. (2017). *Structural vector autoregressive analysis*. Cambridge University Press.  
5. Mertens, K., & Ravn, M. O. (2013). The dynamic effects of personal and corporate income tax changes in the United States. *American Economic Review*, *103*(4), 1212–1247.  
6. Plagborg-Møller, M., & Wolf, C. K. (2021). Local projections and VARs estimate the same impulse responses. *Econometrica*, *89*(2), 955–980.  
7. Ramsay, J. O., & Silverman, B. W. (2005). *Functional data analysis* (2nd ed.). Springer.  
8. Ramsay, J. O., Hooker, G., & Graves, S. (2009). *Functional data analysis with R and MATLAB*. Springer.  
9. Stock, J. H., & Watson, M. W. (2018). Identification and estimation of dynamic causal effects in macroeconomics using external instruments. *The Economic Journal*, *128*(610), 917–948.  

## Citation

If this code contributes to your work, please cite:  
```bibtex
@article{huang2023functional,
  title = {Functional VAR},
  author = {Huang, Jiaming},
  year = {2023},
}
```