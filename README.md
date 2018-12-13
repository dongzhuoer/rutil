# biozhuoer
[![Build Status](https://travis-ci.com/dongzhuoer/biozhuoer.svg?branch=master)](https://travis-ci.com/dongzhuoer/biozhuoer)

## Overview

Some utility functions in bioinformatics, still in developing, use at your own risks.

## Installation

```r
if (!('remotes' %in% .packages(T))) install.packages('remotes');
remotes::install_github('dongzhuoer/biozhuoer');
```

## Usage

refer to `vignette('biozhuoer')`.

## to do

1. add more test for read-fasta (unalign)

## License

[PAUP*](http://phylosolutions.com/paup-test/)
[Aliscore](https://www.zfmk.de/en/research/research-centres-and-groups/aliscore)
[Alicut](https://www.zfmk.de/en/research/research-centres-and-groups/utilities)


## `inst/`

```bash
aria2c -d inst/exec https://raw.githubusercontent.com/PatrickKueck/AliCUT/master/ALICUT_V2.31.pl
```
