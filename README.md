# GeneralizedSasakiNakamura.jl

![license](https://img.shields.io/github/license/ricokaloklo/GeneralizedSasakiNakamura.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/GeneralizedSasakiNakamura.jl.svg)](https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl/releases)
[![Documentation](https://img.shields.io/badge/Documentation-ready)](http://ricokaloklo.github.io/GeneralizedSasakiNakamura.jl)

## Highlights
### Performant frequency-domain Teukolsky solver
Takes on average only a few tens of milliseconds:
<p align="center">
  <img width="50%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248602615-eb611e3f-00b9-49c3-905c-0440b0607683.gif">
</p>

### Solutions that are accurate everywhere
Numerical solutions are *smoothly stitched* to analytical ansatzes near the horizon and infinity at user-specified locations `rsin` and `rsout` respectively:

<p align="center">
  <img width="50%" src="https://github-production-user-asset-6210df.s3.amazonaws.com/55488840/248601145-8ce9b236-d197-4321-83b9-fcb9692a51e8.gif">
</p>

## License
The package is licensed under the MIT License.
