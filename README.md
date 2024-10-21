# Euler inversion: Locating sources of potential-field data through inversion of Euler’s homogeneity equation

by
[Leonardo Uieda](https://leouieda.com),
[Gelson Ferreira Souza-Junior](https://github.com/souza-junior),
[India Uppal](https://github.com/indiauppal),
[Vanderlei Coelho Oliveira Jr.](https://www.pinga-lab.org/people/oliveira-jr.html)

This repository contains the data and source code used to produce the results
presented in:

> Reference of the paper and/or preprint.

|  | Info |
|-:|:-----|
| Version of record | https://doi.org/JOURNAL_DOI |
| Open-access version on EarthArXiv | https://doi.org/EARTHARXIV_DOI |
| Archive of this repository | https://doi.org/10.6084/m9.figshare.26384140 |
| Reproducing our results | [`REPRODUCING.md`](REPRODUCING.md) |

## About

A little bit about this paper, how it came about, and what are the main
contributions. Also include a summary figure or graphical abstract for this
paper.

## Abstract

Earth scientists are able to determine the depth of certain rocks in below ground by measuring the small disturbances that they cause to Earth's gravity and magnetic fields.
A popular method for doing this is *Euler deconvolution*, which is widely available in geoscience software and is fast to run on an average computer. 
Unfortunately, Euler deconvolution has some shortcomings: 1) the approximate shape of the rocks must be known, for example a sphere or a wide flat slab, which is represented by something called the \textit{structural index} 2) the depth of the rocks is not well estimated when there is noise in our data, which is a common occurrence.
We propose a new method that uses more adequate (and complex) mathematics to fix some of the shortcomings of Euler deconvolution.
We call our new method **Euler inversion**.
It is less sensitive to noise in the data and is also able to determine the approximate shape of the source (the structural index).
Our new method is able to replace Euler deconvolution on an Earth scientists toolbox because it is also fast to execute on an average computer.

## License

All Python source code (including `.py` and `.ipynb` files) is made available
under the MIT license. You can freely use and modify the code, without
warranty, so long as you provide attribution to the authors. See
`LICENSE-MIT.txt` for the full license text.

The manuscript text (including all LaTeX files), figures, and data/models
produced as part of this research are available under the [Creative Commons
Attribution 4.0 License (CC-BY)][cc-by]. See `LICENSE-CC-BY.txt` for the full
license text.

[cc-by]: https://creativecommons.org/licenses/by/4.0/
