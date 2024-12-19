# Euler inversion: Locating sources of potential-field data through inversion of Euler’s homogeneity equation

by
[Leonardo Uieda](https://leouieda.com),
[Gelson Ferreira Souza-Junior](https://github.com/souza-junior),
[India Uppal](https://github.com/indiauppal),
[Vanderlei Coelho Oliveira Jr.](https://www.pinga-lab.org/people/oliveira-jr.html)

This repository contains the data and source code used to produce the results
presented in:

> Uieda, L., Souza-Junior, G. F., Uppal, I., Oliveira Jr., V. C. (2024). Euler
> inversion: Locating sources of potential-field data through inversion of
> Euler’s homogeneity equation. EarthArXiv.
> doi:[10.31223/X5T41M](https://doi.org/10.31223/X5T41M).

|  | Info |
|-:|:-----|
| Version of record | TBD |
| Open-access version on EarthArXiv | https://doi.org/10.31223/X5T41M |
| Archive of this repository | https://doi.org/10.6084/m9.figshare.26384140 |
| Reproducing our results | [`REPRODUCING.md`](REPRODUCING.md) |

## About

The main idea for this paper came about during an potential-field methods class
which Leo took in 2012 with his then PhD supervisor [Prof. Valéria C. F.
Barbosa](https://www.pinga-lab.org/people/barbosa.html).
While learning about the Euler deconvolution method, which is a speciality of
Valéria, Leo connected it with the geodetic network adjustment theory he had
been taught by [Prof. Spiros
Pagiatakis](https://www.yorku.ca/spiros/spiros.html) during an exchange program
at York University, Canada, in 2008.
An initial prototype was developed in 2012 but there were still some rough
edges and the project was shelved to make way for other more urgent projects at
the time.
Leo returned to this every few years, making slow progress, and involving
Vanderlei in the planning and discussion of the theory.
In 2024, co-authors Gelson, India, and Vanderlei joined Leo for a sprint to
finish the method and produce this paper.

## Abstract

Earth scientists can estimate the depth of certain rocks beneath Earth's
surface by measuring the small disturbances that they cause in the Earth's
gravity and magnetic fields. A popular method for this is **Euler
deconvolution**, which is widely available in geoscience software and can be
run quickly on a standard computer. Unfortunately, Euler deconvolution has some
shortcomings: 1) the approximate shape of the rocks must be known, for example,
a sphere or a wide flat slab, represented by the **structural index** 2) the
depth of the rocks is not well estimated when there is noise in our data, which
is a common occurrence. We propose a new method, **Euler inversion**, which
fixes some of the shortcomings of Euler deconvolution by using more adequate
(and complex) mathematics. Our method is less sensitive to noise in the data
and is also able to determine the approximate shape of the source (the
structural index). Euler inversion is also fast to execute on a standard
computer, making it a practical alternative to Euler deconvolution on an Earth
scientists toolbox.

<figure>
  <img src="https://github.com/compgeolab/euler-inversion/raw/main/paper/figures/real-data-application.png">
  <figcaption><strong>Figure:</strong> Results of applying Euler inversion with a window size of 12 000 m and a window step of 2400 m to the
aeromagnetic data from Rio de Janeiro, Brazil. Estimated source locations and structural indices obtained from Euler
inversion are shown as triangles (𝜂 = 1), squares (𝜂 = 2), and circles (𝜂 = 3). The colour of each symbol represents
the estimated depth below the surface of the Earth (topography). Also shown are the total-field anomaly flight-line
data, the contours of the post-collisional magmatism and alkaline intrusions (solid black lines) and dykes (dashed
lines). The purple squares highlight the A, B, C, and D anomalies that are discussed in the text.</figcaption>
</figure>

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
