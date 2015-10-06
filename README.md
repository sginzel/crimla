# CRIMLA - CRISPR Modified Landscape Analysis

## Setup

### Requirements
1. Ruby(>= 1.9.2)
2. Git
3. Bundler (>= 1.0.15)
4. Rubygems (>= 1.8.6)

### Installation
```
git clone https://github.com/sginzel/crimla
cd crimla
bundle install
gem build crimla.gemspec
gem install crimla
```

## Usage
```
crimla [options]

  -b, --bam=<s+>                             BAM file(s) to process
  -r, --ref=<s>                              reference file in fasta format
  -c, --coord=<s>                            coordinates to find amplicons
  -t, --tail=<i>                             number of bp to add to coord on both sides (default: 0)
  -e, --exons=<s+>                           Coordinates of the targeted exon (start-stop) (default: )
  -a, --atg=<s>                              Position of the ATG motif in the exon
  -p, --paired, --no-paired                  reads are paired? (default: true)
  -s, --snvs, --no-snvs                      Include SNVs in motifs? (default: true)
  -k, --skip=<i>                             number of reads to skip (default: 0)
  -o, --stop-after=<i>                       Stop processing N reads (default: 0)
  -q, --qnames=<s+>                          List or file of read names to process (default: )
  -m, --min-len=<i>                          minimum length of read (after clipping) (default: 60)
  -i, --mapping-quality=<i>                  include hard clipped bases (default: 50)
  -v, --overlap-only, --no-overlap-only      only consider motifs that are supported by both read pair members. (Default: true)
  -n, --min-overlap=<i>                      Only consider read pairs that overlap at least this many bases. Can be negative to indicate
                                             maximum allowed distance between read pairs. (Default: 0)
  -u, --min-base-quality=<i>                 minimum base quality to consider SNVs (default: 30)
  -f, --min-motif-count=<i>                  minimum count of motif after parsing (default: 10)
  -y, --purity-ploidy=<i>                    Number of expected chromosomes used to calculate the purity index (default: 2)
  -l, --tolerated-indel-len=<i>              maximum length of inframe motifs to still be predicted harmless. Do not use this to make your
                                             knock-out look more successful. Use it ONLY if you know what you are doing. (Default: -1)
  -g, --coverage, --no-coverage              calculate coverages per motif (default: true)
  --coverage-global, --no-coverage-global    Show global coverage (default: true)
  --coverage-motif                           Show coverage for every motif (recommended for single-end reads)
  --top=<i>                                  show top N motifs (default: 50)
  --csv=<s>                                  filename of the CSV file to produce
  --svg=<s>                                  filename of the SVG file to produce
  --version                                  Print version and exit
  -h, --help                                 Show this message
```

## Example
```
crimla \
 --bam ITK1.bam \
 --ref human.fasta \
 --svg ITK1.svg \
 --csv ITK1.csv \
 --min-len 80 \
 --mapping-quality 50 \
 --top 25 \
 --tail 25 \
 --coord 5:156,608,004-156,608,126 \
 --exons 5:156,607,837-156,608,126 \
 --atg 156,607,989 \
 --paired \
 --overlap-only \
 --coverage \
 --coverage-global \
 --no-coverage-motif \
 --snvs
```
