# a-continuous-time-stochastic-Boolean-model-of-the-budding-yeast-cell-cycle

## 'model.R'
- requires 'pars.R' (parameter values of the model).
- produces a variable named 'state' which is a dataframe containing a sequence of 14 recurring states shown in Fig. 1b and Table 3.

## 'simulate Fig2.R'
- requires 'pars.R' (parameter values of the model).
- produces a PDF file of Fig 2 in the manscript.

## 'simulate Fig3 & 4.R'
- requires 'pars.R' (parameter values of the model).
- requires 'stat_exp.csv' (experimental data from <a href="https://www.nature.com/articles/nature06072" target="_blank">Di Talia et. al's 2007 paper</a>)
- produces PDF files of Figs 3 and 4 in the manuscript.
