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
1. move `read_sam()`, `read_bed()`, etc from paristools to here.