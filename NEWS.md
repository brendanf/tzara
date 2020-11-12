# tzara 0.0.11

* fixed a bug in `cluster_consensus()` with character inputs when names are
  not given.
* `cluster_consensus()` now gives a result for a single (non-`NA`) unique
  sequence, even if there are less than 3 copies.

# tzara 0.0.10

* `extract_region()` auto-detects the format of the output file, if given,
  from the name.
* Fixed a bug where `extract_region()` was not saving quality information
  from the input to an output "fastq" file when the input was a "fastq" file or
  `QualityScaledXStringSet`.
* Added argument `"dir"` to `dadamap()`, to specify the location of input files,
  if `derepFastq()` was originally called on a directory of files.
* fixed a bug with `"seq_id"` and `"filename"` not being passed correctly from
  `dadamap.list()` to `dadamap.derep`.
* `has_alphabet()` now works even if some of the input sequences are `NA`.

# tzara 0.0.9

* Added arguments `"seq_id"` and `"filename"` to `dadamap()`, to make it more
  user-friendly.  Use of `add_dada_names()` is now unnecessary.
  
# tzara 0.0.8

* Added a `NEWS.md` file to track changes to the package.
* **Breaking Change**: Uniformly use the column name "seq_id" rather than "seq",
  "seq.id", or "seq_name" for unique sequence identifiers.
* Add support for `DNAStringSet` and `QualityScaledDNAStringSet` in 
  `extract_region()`.
