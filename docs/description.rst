Description of the methods
==========================

Viral genome analysis
---------------------

Viral genome assembly
~~~~~~~~~~~~~~~~~~~~~

The filtered and trimmed reads are subsampled to at most 100,000 pairs.
*de novo* assemby is performed using Trinity_. SPAdes_ is also offered as
an alternative *de novo* assembler.
Reference-assisted assembly improvements follow (contig scaffolding, orienting, etc.)
with MUMMER_ and MUSCLE_ or MAFFT_. Gap2Seq_ is used to seal gaps between scaffolded *de novo* contigs with sequencing reads.

Each sample's reads are aligned to its *de novo* assembly using Novoalign_
and any remaining duplicates were removed using Picard_ MarkDuplicates.
Variant positions in each assembly were identified using GATK_ IndelRealigner and
UnifiedGenotyper on the read alignments. The assembly was refined to represent the
major allele at each variant site, and any positions supported by fewer than three
reads were changed to N.

This align-call-refine cycle is iterated twice, to minimize reference bias in the assembly.
 
.. _Trinity: http://trinityrnaseq.github.io/
.. _SPAdes: http://bioinf.spbau.ru/en/spades
.. _MUMMER: https://mummer4.github.io/
.. _MUSCLE: https://www.drive5.com/muscle/
.. _MAFFT: http://mafft.cbrc.jp/alignment/software/
.. _Gap2Seq: https://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/
.. _Novoalign: http://www.novocraft.com/products/novoalign/
.. _Picard: http://broadinstitute.github.io/picard
.. _GATK: https://www.broadinstitute.org/gatk/

