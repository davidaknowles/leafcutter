The leafcutter R package is a subdir of the leafcutter repo which makes automatically building the conda recipe difficult. A workaround is to make a new github repo with just the R package: https://github.com/davidaknowles/lc_conda (note this must have at least one `tag`) Then run
```conda skeleton cran https://github.com/davidaknowles/lc_conda```

Additional steps: 
* change the `git_url` to point to the main repo and `git_tag` to the latest tag `v0.2.1`
* edit `build.sh` to include `--preclean` so that `stan` regenerates the C++ code, and to operate in the `leafcutter` subdir of the source dir
* add 'r-roxygen2' as a build dependency in `meta.yaml`, as this is required by the preclean script

Not the Windows `bld.bat` is included here in case anyone wants to try building on Windows, but it's untested. 

To build the conda package: `conda build .` 
To install locally: `conda install -c local r-leafcutter`

My full conda setup on linux-64 was:
```
wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -bfp /scr/miniconda2/
. /scr/miniconda2/bin/activate
conda install -y conda-build
conda update -y --all
conda install -y -c r r-essentials=1.6.0
export CONDA_R=3.4.1

conda skeleton cran intervals
conda-build r-intervals
conda install -y -c local r-intervals

conda skeleton cran optparse
conda-build r-optparse
conda install -y -c local r-optparse

conda skeleton cran shinycssloaders
conda-build r-shinycssloaders
conda install -y -c local r-shinycssloaders 

conda skeleton cran https://github.com/davidaknowles/lc_conda
```
(manual edits as described above)

```
conda install -y -c local r-leafcutter
```

