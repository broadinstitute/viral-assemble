# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

viral-assemble is a set of scripts and tools for the assembly of viral genomes from NGS data, primarily used for Lassa and Ebola virus analysis. This is a docker-centric Python project built on top of viral-core, with a modular architecture for viral genome assembly workflows.

## Development Commands

### Testing

Run all unit tests:
```bash
pytest -rsxX -n auto test/unit
```

Run specific test file:
```bash
pytest test/unit/test_assembly.py
```

Run integration tests:
```bash
pytest test/unit/test_assembly_integration.py
```

Run with coverage:
```bash
pytest --cov
```

### Docker Development Workflow

The development paradigm is intentionally docker-centric. To work on code changes:

1. Mount local checkout into viral-core container:
```bash
docker run -it --rm -v $(pwd):/opt/viral-ngs/viral-assemble quay.io/broadinstitute/viral-core
```

2. If changing conda dependencies, update them inside container:
```bash
/opt/viral-ngs/source/docker/install-conda-dependencies.sh /opt/viral-ngs/viral-assemble/requirements-conda.txt
```

3. Or install full dev layer from viral-core:
```bash
/opt/viral-ngs/viral-assemble/docker/install-dev-layer.sh
```

4. Test interactively within container:
```bash
cd /opt/viral-ngs/viral-assemble
pytest -rsxX -n auto test/unit
```

### Docker Build

Build docker image:
```bash
docker build -t viral-assemble .
```

The Dockerfile layers viral-assemble on top of viral-core:2.4.2, installing conda dependencies first, then copying source code.

## Architecture

### Main Entry Point

`assembly.py` - Main command-line interface with multiple subcommands registered via the `__commands__` list. Uses argparse for CLI and util.cmd for command registration.

### Core Assembly Commands

Key subcommands available via `assembly.py <command>`:
- `assemble_spades` - De novo assembly using SPAdes
- `refine_assembly` - Iterative refinement using reference-guided assembly
- `order_and_orient` - Order and orient contigs against reference
- `impute_from_reference` - Fill gaps using reference sequence
- `gapfill_gap2seq` - Close gaps between contigs using Gap2Seq
- `cluster_references_ani` - Cluster reference genomes by ANI
- `skani_contigs_to_refs` - Find closest references for assembled contigs
- `vcf_to_fasta` - Convert VCF to consensus FASTA
- `trim_rmdup_subsamp` - Read preprocessing pipeline
- `filter_short_seqs` - Remove short sequences from assembly
- `modify_contig` - Modify contig sequences programmatically

### Module Structure

- `assemble/` - Tool wrapper modules for assembly algorithms
  - `spades.py` - SPAdes assembler wrapper
  - `mafft.py` - MAFFT multiple sequence aligner
  - `mummer.py` - MUMmer alignment and variant calling
  - `muscle.py` - MUSCLE aligner wrapper
  - `gap2seq.py` - Gap2Seq gap filling
  - `skani.py` - skani ANI-based clustering and comparison
  - `vcf.py` - VCF parsing and variant calling utilities

- `test/` - pytest-based test suite
  - `test/unit/test_assembly.py` - Unit tests for assembly functions
  - `test/unit/test_assembly_integration.py` - Integration tests for full workflows
  - `conftest.py` - pytest fixtures and configuration

### Dependencies from viral-core

viral-assemble imports core utilities from viral-core:
- `util.cmd` - Command-line parsing and command registration
- `util.file` - File handling utilities
- `util.misc` - Miscellaneous utilities
- `read_utils` - Read processing utilities
- `tools.*` - Tool wrappers (picard, samtools, gatk, novoalign, trimmomatic, minimap2)

### Conda Dependencies

Tool dependencies are specified in `requirements-conda.txt` and installed via conda:
- gap2seq, mafft, mummer4, muscle, spades - Assembly/alignment tools
- skani, fastani, sourmash - ANI comparison and sketching
- rasusa - Downsampling
- perl-bio-easel, sequip - Perl bioinformatics tools

## Testing Requirements

- pytest is used with parallelized execution (`-n auto`)
- Tests use fixtures from `test.stubs` and `conftest.py`
- Test input files are in `test/input/<TestClassName>/`
- Access test inputs via `util.file.get_test_input_path(self)`
- **New tests should add no more than ~20-30 seconds to testing time**
- **Tests taking longer must be marked with `@pytest.mark.slow`**
- Run slow tests with `pytest --runslow`
- **New functionality must include unit tests covering basic use cases and confirming successful execution of underlying binaries**

## CI/CD

GitHub Actions workflow (`.github/workflows/build.yml`) runs on push/PR:
- Docker image build and push to quay.io/broadinstitute/viral-assemble
  - Master branch: tagged as `latest` and with version number
  - Non-master branches: tagged with branch name (ephemeral, deleted within two weeks)
- Unit and integration tests with pytest
- Documentation build validation
- Conda package build and deploy

## Key Design Patterns

### Command Registration
Commands are registered by appending `(command_name, parser_function)` tuples to `__commands__`. Each command has:
- A parser function (`parser_<command_name>`) that creates argparse parser
- A main function that implements the logic
- Connection via `util.cmd.attach_main(parser, main_function, split_args=True)`

### Assembly Pipeline Flow
Typical assembly workflow:
1. `trim_rmdup_subsamp` - Clean and subsample reads
2. `assemble_spades` - De novo assembly
3. `order_and_orient` - Scaffold against reference
4. `impute_from_reference` - Fill gaps from reference
5. `refine_assembly` - Iteratively improve (combines multiple steps)

### Error Handling
- `DenovoAssemblyError` - Raised for assembly failures or insufficient data
- `IncompleteAssemblyError` - Raised when assembly quality thresholds not met
- `PoorAssemblyError` - Raised for assemblies failing quality criteria
