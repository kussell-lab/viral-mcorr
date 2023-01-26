# viral-mcorr
viral-mcorr is a method for inferring recombination rates from 
large-scale sequencing data in (+)ssRNA viruses using correlation profiles of synonymous substitutions.

The viral-mcorr method is described in the following paper:

```bibtex
@article {doi:10.1073/pnas.2206945119,
    author = {Asher Preska Steinberg  and Olin K. Silander  and Edo Kussell },
    title = {Correlated substitutions reveal SARS-like coronaviruses recombine frequently with a diverse set of structured gene pools},
    journal = {Proceedings of the National Academy of Sciences},
    volume = {120},
    number = {5},
    pages = {e2206945119},
    year = {2023},
    doi = {10.1073/pnas.2206945119},
    URL = {https://www.pnas.org/doi/full/10.1073/pnas.2206945119}
}
```
[https://www.pnas.org/doi/full/10.1073/pnas.2206945119](https://www.pnas.org/doi/full/10.1073/pnas.2206945119)

# Requirements
* Install `git` from [https://git-scm.com](https://git-scm.com/);
* Install `go` from [https://golang.org/doc/install](https://golang.org/doc/install);
* Install `python3` from [https://www.python.org/](https://www.python.org/) (we found running issues using the default Python in MacOS);
* Install `pip3` from [https://pip.pypa.io/en/stable/installing/](https://pip.pypa.io/en/stable/installing/).

# Installation
1. For basic usage, install `mcorr-gene-aln`, `mcorrViralGenome`, `mcorrLDGenome` from your terminal:
```sh

go install github.com/kussell-lab/viral-mcorr/cmd/mcorr-gene-aln@latest
go install github.com/kussell-lab/viral-mcorr/cmd/mcorrViralGenome@latest
go install github.com/kussell-lab/viral-mcorr/cmd/mcorrLDGenome@latest


cd $HOME/go/src/github.com/kussell-lab/mcorr/cmd/mcorr-viral-fit
pip install $HOME/go/src/github.com/kussell-lab/mcorr/cmd/mcorr-viral-fit
```
Install `mcorr-viral-fit` by cloning this github repository and then using pip to install the program locally:

```sh
git clone git@github.com:kussell-lab/viral-mcorr.git
pip install ./
```
2. Add `$HOME/go/bin` and `$HOME/.local/bin` to your `$PATH` environment. In Linux, you can do it in your terminal:
```sh
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin
```

In MacOS, you can do it as follows:
```sh
export PATH=$PATH:$HOME/go/bin:$HOME/Library/Python/3.6/bin
```

We have tested installation in MacOS Monterey (w/ an M1 chip), using Python 3 and Go 1.15 and 1.16.

## Basic usage for inferring recombination parameters
The inference of recombination parameters requires two steps:

1. Calculate _Correlation Profile_

    1. For multi-fasta alignments of single genes or whole genomes in which
       there is a single CDS region, use `mcorr-gene-aln` :
       ```sh
       mcorr-gene-aln <input MFA file> <output prefix>
       ```
       
    2. To calculate correlation profiles across the CDS region of whole-genome alignments (multiple gene alignments), use `mcorrViralGenome`:

       ```sh
       mcorrViralGenome <input XMFA file> <output prefix>
       ```
       The flag `--mate-aln` allows for inclusion of a second XMFA file of viral genomes. 
       The flag `--between-clades` can be used when you have two XMFA files to calculate correlation profiles exclusively across
       sequence pairs in which neither sequence is from the same XMFA file.
       
        The XMFA files should contain only *coding* sequences and should not include any redundant CDS regions 
       (i.e., CDS regions which code for a subregion of another CDS region should be removed from the XMFA). Gapped regions should be denoted by dashes or Ns. 
       The description of XMFA file can be found in [http://darlinglab.org/mauve/user-guide/files.html](http://darlinglab.org/mauve/user-guide/files.html). We provide two useful pipelines to generate whole-genome alignments:
        * from multiple assemblies: [https://github.com/kussell-lab/AssemblyAlignmentGenerator](https://github.com/kussell-lab/AssemblyAlignmentGenerator);
        * from raw reads: [https://github.com/kussell-lab/ReferenceAlignmentGenerator](https://github.com/kussell-lab/ReferenceAlignmentGenerator)
    

   All programs will produce two files:
    * a .csv file stores the calculated Correlation Profile, which will be used for fitting in the next step;
    * a .json file stores the (intermediate) Correlation Profile for each gene.

2. Fit the Correlation Profile using `mcorr-viral-fit`:
    1. For fitting correlation profiles as described in our paper [link will go here] use `mcorr-viral-fit`:

          ```sh
          mcorr-viral-fit <.csv file> <output_prefix>
          ```

       This will produce several files:

        * `<output_prefix>_template-switch_best_fit.svg` and `<output_prefix>_zero-recombo_best_fit.svg` show the plots of the Correlation Profile, fitting, and residuals for the template-switching recombination model and for the zero recombination case;
        * `<output_prefix>_comparemodels.csv` shows the table of fitted parameters for all recombination models (template-switching, fragment-incorporation, and zero-recombination) and AIC values;
        * `<output_prefix>_template-switch_residuals.csv` and `<output_prefix>_zero-recombo_residuals.csv` includes residuals for the model with template-switching and the zero-recombination case
        * `<output_prefix>_template-switch_fit_results.csv` shows fit results for data and bootstrap replicates to template-switching model (if correlations were analyzed w/ `mcorr-gene-aln`)
        * `<output_prefix>_template-switch_fit_report.txt` shows fit results and bootstrap CIs if correlations were analyzed w/ `mcorr-gene-aln`

## Basic usage for measuring correlation coefficients for sites across the genome or genes
To measure correlations at individual codons across the genome, you can use `mcorrLDGenome` as 
described in our paper [link will go here]:


          mcorrLDGenome <input XMFA file> <output prefix>

XMFA files must be formatted in the same way as described for mcorrViralGenome, above. Alternatively, multi-fasta alignments of 
single CDS regions.

# Examples

1. [How to create alignments of viral genomes for use with viral-mcorr.](https://github.com/kussell-lab/virus_alignment_example)
2. [Example workflows for using mcorr-gene-aln, mcorrViralGenome, and mcorrLDGenome](https://github.com/kussell-lab/viral-mcorr_sl-cov_examples)


    

