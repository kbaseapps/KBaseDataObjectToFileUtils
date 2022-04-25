### Version 1.1.1
__Changes__
- adjust fasta IDs to be under length threshold and return ID map
- change STOP residue from '_' to '*'
- support genetic codes 1, 2, 3, 4, 5, 6, and 11
- support feature-specific genetic code (not yet available)
- don't write protein translations with internal STOP codon

### Version 1.1.0
__Changes__
- added SpeciesTreeToFASTA() method

### Version 1.0.1
__Changes__
- include 'non_coding_features' for nucleotide/non-CDS queries from Genomes >= 9.0
- added more tests for nucleotide mode

### Version 1.0.0
__Changes__
- handle gene annotation functions from Genomes >= 9.0 (Features >= 3.0)
- support FeatureSets with AnnotatedMetagenomeAssembly features
- added Github Actions testing
- made all test data uploaded

### Version 0.0.9
__Changes__
- Added genome_ref_to_obj_name return value

### Version 0.0.8
__Changes__
- Added AnnotatedMetagenomeAssemblyToFASTA() method

### Version 0.0.7
__Changes__
- patched CheckJob() "Bad Status Line" error
